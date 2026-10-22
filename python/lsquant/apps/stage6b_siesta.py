"""Stage-6b gate executor; the CLI only constructs a stage-7 RunFile.

No material data or expected plateau is embedded here. The same loaded
SIESTA matrices feed the Bloch engine, eigenphase reference, and NumPy
oracle. Nothing is certified as a golden by this program.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

from .. import _core
from ..bases import BlochMetricBasis
from ..models.siesta import load_siesta_model
from ..reconstruct import dos_curve, shc
from ..words import compile as compile_word, parse as parse_word
from ..words.grammar import OperatorSlot
from ..oracle.spin_hall import UNITS, spin_berry_curvature, sigma_xy, spin_chern

SHC_WORD = "Tr[JXSZ delta(E-H) VY dG-(E)/dE]"
DOS_WORD = "Tr[delta(E-H)]"


def compute_moments(basis, word, *, moments=1024, trace="random", states=16,
                    seed=12345, window=None, batch=None):
    """Compile the physical word and run the unchanged C++ moment engine.

    For a nonorthogonal Bloch calculation the compiler builds the stage-6
    adjoint bra ``S T_m X A X r`` and ket ``X B T_n r``. bind_eigen_trace
    replaces the two seeds by SC r and C r, respectively, so each sample
    equals the eigenbasis random-phase sample. Rank-2 weights are applied
    once without imposing moment symmetry on a finite stochastic trace.
    """
    eigen = basis.name == "eigen"
    parsed = parse_word(word)
    operators = {slot.operator_key: basis.operator(slot.operator_key)
                 for slot in parsed.slots if isinstance(slot, OperatorSlot)}
    if window is None:
        window = basis.window()
    generator = basis.states(trace, None if trace == "exact" else states, seed)
    nstates = int(generator.number_of_states)
    metric = ({"policy": "orthogonal", "route": "lowdin"} if eigen or basis.route == "anticommutator" else
              {"policy": "nonorthogonal", "route": basis.route,
               "applier": "exact", "metric_object": basis.metric})
    compiled = compile_word(parsed, basis={"name": basis.name, "operators": operators,
                                          "metric": metric},
                            params={"expansion": {"moments": moments, "states": nstates,
                                                   "window": window}, "metric": metric})
    if compiled.match.entry.shape_id not in ("dos_identity", "dos_operator", "kubo_bastin"):
        raise ValueError("stage6b supports spectral and Bastin words only")
    plan, engine_ops = compiled.plan, compiled.operators
    if eigen:
        h = _core.RescaledOperator(basis.hamiltonian(), window)
    else:
        h = basis.rescaled_hamiltonian(window)
        plan, engine_ops = basis.bind_eigen_trace(plan, engine_ops)
    engine = _core.MomentEngine(h, engine_ops, generator, moments if batch is None else batch)
    mu = engine.run(plan)
    mu.window, mu.system_size = window, basis.dimension
    mu.number_of_states, mu.state_label, mu.word = nstates, generator.label, compiled.plan.word
    if mu.rank == 2:
        _core.ApplyKpmWeights(mu, nstates, False)
    else:
        mu.kpm_weighted = True
    return mu, compiled


def _write_json(path, data):
    """Write finite JSON, refusing NaN/Infinity instead of hiding an invalid oracle."""
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def execute_siesta_runfile(rf, out_dir):
    """Execute a validated Siesta RunFile and persist its resolved provenance.

    Normalisation: ``SHC=2*pi*Bastin/(N_k*A_cell)``; DOS is per cell.
    The eigen reference uses identical phase seed/count and the same exact
    S^-1/2 factorisation. The production and reference moment differences
    are reported without declaring a tolerance verdict or a gate result.
    """
    if rf.basis["name"] not in ("bloch", "eigen"):
        raise ValueError("SIESTA stage6b requires basis.name bloch or eigen")
    if rf.metric["policy"] != "nonorthogonal" or rf.metric["applier"] != "exact":
        raise ValueError("SIESTA stage6b requires nonorthogonal metric with exact per-block applier")
    if rf.basis["disorder"]["W"] != 0:
        raise ValueError("nonorthogonal disorder belongs to stage6c")
    if rf.expansion["trace"] != "random":
        raise ValueError("material comparisons require eigenbasis random-phase trace")
    params = rf.model["params"]
    mode = params.get("velocity_mode", "berry_connection")
    model = load_siesta_model(rf.model["dir"], rf.model["label"], dipoles=mode == "covariant",
                              quadrature_order=int(params.get("quadrature_order", 24)),
                              quadrature_tolerance=float(params.get("quadrature_tolerance", 1e-5)),
                              quadrature_max_order=int(params.get("quadrature_max_order", 96)))
    bloch = BlochMetricBasis(model, rf.basis["supercell"], velocity_mode=mode,
                             route=rf.metric["route"], seed=rf.expansion["seed"])
    basis = bloch if rf.basis["name"] == "bloch" else bloch.eigenbasis()
    window_spec = rf.expansion["window"]
    if window_spec == "gershgorin":
        raise ValueError("SIESTA windows must come from the generalised block eigenvalues")
    window = basis.window() if window_spec == "eigen" else _core.SpectralWindow(*window_spec)
    if window.lower > bloch.eigenvalues.min() or window.upper < bloch.eigenvalues.max():
        raise ValueError("explicit window does not enclose the generalised spectrum")
    area = float(np.linalg.norm(np.cross(model.uc[0], model.uc[1])))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    resolved = rf.to_resolved_dict()
    resolved["output"]["dir"] = str(out)
    resolved["reconstruction"].update(num_cells=bloch.num_cells, cell_area=area)
    # Keep resolved JSON valid as a RunFile: extra provenance lives separately.
    resolved["expansion"]["window"] = [float(window.lower), float(window.upper)]
    _write_json(out/"run.resolved.json", resolved)
    _write_json(out/"basis.json", dict(model.metadata, metric_residuals=bloch.metric_residuals,
                                      window=[float(window.lower), float(window.upper)],
                                      velocity_mode=mode, route=bloch.route,
                                      metric_approximation=bloch.route == "anticommutator",
                                      trace_basis="generalised eigenvectors" if bloch.route != "anticommutator"
                                      else "eigenvectors of the approximate anticommutator Hamiltonian",
                                      units=UNITS, num_cells=bloch.num_cells, cell_area=area))
    moment_files, curve_files, comparisons = [], [], []
    broadening = rf.reconstruction["broadening_meV"]
    if broadening is None:
        raise ValueError("stage6b requires reconstruction.broadening_meV")
    states = rf.expansion["states"] or 16
    supported = {parse_word(SHC_WORD).canonical(), parse_word(DOS_WORD).canonical()}
    if any(parse_word(word).canonical() not in supported for word in rf.words):
        raise ValueError("stage6b run files accept the physical SHC word or identity DOS word")
    allowed_formulas = {"shc", "dos", "dos_curve"}
    if rf.reconstruction["formulas"] and not set(rf.reconstruction["formulas"]) <= allowed_formulas:
        raise ValueError("stage6b reconstruction supports shc and dos")
    if rf.reconstruction["targets"]:
        raise ValueError("stage6b does not accept legacy tutorial reconstruction targets")
    for word in rf.words:
        mu, compiled = compute_moments(basis, word, moments=rf.expansion["moments"], states=states,
                                       seed=rf.expansion["seed"], window=window, batch=rf.expansion["batch"])
        if mu.rank == 2 and parse_word(word).canonical() != parse_word(SHC_WORD).canonical():
            raise ValueError("the shc formula requires the physical JXSZ/VY Bastin word")
        prefix = compiled.reconstruction.legacy_moment_prefix
        size = f"{rf.expansion['moments']}" + (f"x{rf.expansion['moments']}" if mu.rank == 2 else "")
        stem = f"{prefix}{model.seed_name}KPM_M{size}_state{mu.state_label}"
        path = out/(stem+".chebmom")
        _core.WriteChebmom(str(path), mu)
        moment_files.append(str(path))
        if rf.output["legacy_layout"]:
            path = out/(stem+compiled.reconstruction.legacy_moment_ext)
            writer = _core.WriteLegacyChebmom2D if mu.rank == 2 else _core.WriteLegacyChebmom1D
            writer(str(path), mu)
            moment_files.append(str(path))
        # Preserve full precision unless the run file explicitly requests the legacy round trip.
        reconstruction_mu = (_core.ReadChebmom(str(path))
                             if rf.reconstruction["legacy_text_round_trip"] else mu)
        if mu.rank == 2:
            energy, values = shc(reconstruction_mu, broadening, bloch.num_cells, area)
            curve_name, header = "shc", f"energy_eV sigma_z_xy [{UNITS}]"
        else:
            energy, values = dos_curve(reconstruction_mu, broadening)
            values = values/bloch.num_cells
            curve_name, header = "dos", "energy_eV DOS [states/eV/cell]"
        if bloch.route == "anticommutator":
            header += "; approximate anticommutator metric"
        curve_path = out/(curve_name+".dat")
        np.savetxt(curve_path, np.column_stack((energy, values)), fmt="%.17g", header=header)
        curve_files.append(str(curve_path))
        # The independent engine comparison uses a diagonal H and C^+ O C,
        # not the covariant matrices passed through a second metric route.
        reference, _ = compute_moments(bloch.eigenbasis(), word, moments=rf.expansion["moments"],
                                        states=states, seed=rf.expansion["seed"], window=window,
                                        batch=rf.expansion["batch"])
        reference_path = out/("eigen_"+curve_name+".chebmom")
        _core.WriteChebmom(str(reference_path), reference)
        moment_files.append(str(reference_path))
        if mu.rank == 2:
            er, vr = shc(reference, broadening, bloch.num_cells, area)
        else:
            er, vr = dos_curve(reference, broadening)
            vr = vr/bloch.num_cells
        eigen_curve = out/("eigen_"+curve_name+".dat")
        np.savetxt(eigen_curve, np.column_stack((er, vr)), fmt="%.17g", header=header)
        curve_files.append(str(eigen_curve))
        comparisons.append(dict(word=word, phase_seed=rf.expansion["seed"], states=states,
                                max_moment_difference=float(np.max(np.abs(mu.data-reference.data))),
                                max_curve_difference=float(np.max(np.abs(values-vr)))))
    _write_json(out/"eigen_comparison.json", comparisons)
    if params.get("oracle", False):
        curvature = spin_berry_curvature(bloch.physical_eigenvalues, bloch.physical_eigenvectors,
                                         bloch.operator_blocks("JXSZ"), bloch.operator_blocks("VY"))
        energies = np.linspace(float(window.lower), float(window.upper),
                               rf.reconstruction["num_points"] or 1001)
        curve = sigma_xy(bloch.physical_eigenvalues, curvature, energies, area, broadening=broadening/1000)
        _write_json(out/"exact_shc.json", dict(energies=energies.tolist(), values=curve.tolist(), units=UNITS,
                                              trace="eigenbasis random phase (diagonal curvature trace is exact)",
                                              states=states, seed=rf.expansion["seed"],
                                              gaussian_standard_deviation_eV=broadening/1000))
        ef = (rf.reconstruction["fermi_energies"] or [0.0])[0]
        try:
            topology = spin_chern(bloch.physical_eigenvalues, bloch.physical_eigenvectors, bloch.sk_array,
                                  bloch.operator_blocks("SZ"), ef, bloch.supercell,
                                  k_points=bloch.k_points, centres=model.tau)
        except ValueError as error:
            # A metallic input is a scientific result, not an invented integer.
            topology = {"defined": False, "reason": str(error)}
        topology.update(fermi_energy_eV=ef,
                        sigma_xy=float(sigma_xy(bloch.physical_eigenvalues, curvature, [ef], area,
                                                 broadening=broadening/1000)[0]), units=UNITS)
        _write_json(out/"spin_chern.json", topology)
    return dict(output_dir=str(out), resolved_runfile=str(out/"run.resolved.json"),
                moment_files=moment_files, curve_files=curve_files,
                written_files=[str(p) for p in sorted(out.iterdir()) if p.is_file()])


def main(argv=None):
    """Translate command-line choices into the same RunFile used by lsquant run."""
    from ..runfile import execute_runfile
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--siesta", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--supercell", nargs=3, type=int, default=[48, 48, 1])
    parser.add_argument("--velocity-mode", choices=["bare", "berry_connection", "covariant"], default="berry_connection")
    parser.add_argument("--route", choices=["lowdin", "generalised", "anticommutator"], default="lowdin")
    parser.add_argument("--word", choices=["shc", "dos"], default="shc")
    parser.add_argument("--moments", type=int, default=1024)
    parser.add_argument("--trace", choices=["random"], default="random")
    parser.add_argument("--states", type=int, default=16)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--broadening-meV", type=float, default=10.)
    parser.add_argument("--out", required=True)
    parser.add_argument("--oracle", action="store_true")
    parser.add_argument("--quadrature-order", type=int, default=24)
    parser.add_argument("--quadrature-max-order", type=int, default=96)
    parser.add_argument("--quadrature-tolerance", type=float, default=1e-5)
    args = parser.parse_args(argv)
    raw = dict(model=dict(kind="siesta", dir=str(Path(args.siesta).expanduser().resolve()), label=args.label,
                          params=dict(velocity_mode=args.velocity_mode, oracle=args.oracle,
                                      quadrature_order=args.quadrature_order,
                                      quadrature_max_order=args.quadrature_max_order,
                                      quadrature_tolerance=args.quadrature_tolerance)),
               basis=dict(name="bloch", supercell=args.supercell),
               metric=dict(policy="nonorthogonal", route=args.route, applier="exact"),
               word=SHC_WORD if args.word == "shc" else DOS_WORD,
               expansion=dict(moments=args.moments, trace=args.trace, states=args.states,
                              seed=args.seed, window="eigen"),
               reconstruction=dict(formulas=[args.word], broadening_meV=args.broadening_meV),
               output=dict(dir=args.out, legacy_layout=True, write_chebmom=True))
    print(json.dumps(execute_runfile(raw), indent=2))


if __name__ == "__main__":
    main()
