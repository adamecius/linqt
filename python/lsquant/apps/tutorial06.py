"""Reproduce Tutorial 06 (Elliott-Yafet and Dyakonov-Perel spin relaxation on a 2D square lattice) through the Python API.

Imports ``examples/06_spin_relaxation_square/make_ey.py`` and
``examples/06_spin_relaxation_square/make_rashba.py``, evaluates the projected
spin time-evolution words ``Tr[P:PZ U(t)^+ SZ U(t) delta(E-H) P:PZ]``,
``Tr[P:PZ U(t)^+ 1 U(t) delta(E-H) P:PZ]``, ``Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]``,
and ``Tr[P:PX U(t)^+ 1 U(t) delta(E-H) P:PX]`` through :func:`~lsquant.words.compile`
and :class:`~lsquant._core.MomentEngine`, reconstructs the time-correlation
curves at ``E_F = 2.0`` eV (``140`` meV Jackson broadening), and writes the
top-level ``.dat`` files, ``ey_ratelaw.dat``, and ``captured/0000..0254``.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import importlib.util
import os
import shutil
import sys
import tempfile
from types import ModuleType

import numpy as np
import numpy.typing as npt

from .. import _core
from ..moments import run_plan
from ..reconstruct import time_correlation_curve
from ..words import compile as compile_word


def _load_tutorial06_generators() -> tuple[ModuleType, ModuleType]:
    """Import ``make_ey.py`` and ``make_rashba.py`` from ``examples/06_spin_relaxation_square``."""
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..")
    )
    tut_dir = os.path.join(repo_root, "examples", "06_spin_relaxation_square")
    ey_path = os.path.join(tut_dir, "make_ey.py")
    rashba_path = os.path.join(tut_dir, "make_rashba.py")
    if not os.path.isfile(ey_path):
        raise FileNotFoundError(f"tutorial06: make_ey.py not found at {ey_path}")
    if not os.path.isfile(rashba_path):
        raise FileNotFoundError(
            f"tutorial06: make_rashba.py not found at {rashba_path}"
        )

    spec_ey = importlib.util.spec_from_file_location("make_ey", ey_path)
    if spec_ey is None or spec_ey.loader is None:
        raise RuntimeError(f"tutorial06: could not load spec for {ey_path}")
    mod_ey = importlib.util.module_from_spec(spec_ey)
    sys.modules["make_ey"] = mod_ey
    spec_ey.loader.exec_module(mod_ey)

    spec_ra = importlib.util.spec_from_file_location("make_rashba", rashba_path)
    if spec_ra is None or spec_ra.loader is None:
        raise RuntimeError(f"tutorial06: could not load spec for {rashba_path}")
    mod_ra = importlib.util.module_from_spec(spec_ra)
    spec_ra.loader.exec_module(mod_ra)
    return mod_ey, mod_ra


def _read_bounds(bounds_path: str) -> _core.SpectralWindow:
    """Read ``lo hi`` from a legacy ``BOUNDS`` file."""
    with open(bounds_path, encoding="utf-8") as fh:
        tokens = fh.read().split()
    return _core.SpectralWindow(float(tokens[0]), float(tokens[1]))


def _fit_rate(
    t: npt.NDArray[np.float64],
    s: npt.NDArray[np.float64],
    lo: float = 0.2,
    hi: float = 0.85,
) -> float:
    """Fit exponential spin relaxation time ``T_1`` over ``lo < s < hi``."""
    mask = (s > lo) & (s < hi)
    if int(np.sum(mask)) < 3:
        return float("nan")
    slope = float(np.polyfit(t[mask], np.log(s[mask]), 1)[0])
    return -1.0 / slope


def _eval_and_reconstruct_pair(
    stage_dir: str,
    label: str,
    op_name: str,
    proj_name: str,
    num_moments: int,
    ntime: int,
    tmax: float,
    states: _core.StateGenerator,
    state_label: str,
    broadening_meV: float,
    ef: float,
    bessel: str,
    legacy_text_round_trip: bool,
    out_dir: str,
    captured_dir: str,
    cap_idx: int,
    *,
    keep_top_level: bool = True,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64], list[str], int]:
    """Run numerator (``op_name``) and denominator (``1``) projected words and reconstruct at ``ef``."""
    op_dir = os.path.join(stage_dir, "operators")
    window = _read_bounds(os.path.join(stage_dir, "BOUNDS"))
    ham = _core.ReadCsrText(os.path.join(op_dir, f"{label}.HAM.CSR"))
    dim = int(ham.dimension)

    proj_mat = _core.ReadCsrText(os.path.join(op_dir, f"{label}.{proj_name}.CSR"))
    proj_op = _core.SparseMatrixOperator(proj_mat)
    op_mat = _core.ReadCsrText(os.path.join(op_dir, f"{label}.{op_name}.CSR"))
    op_lin = _core.SparseMatrixOperator(op_mat)
    id_op = _core.IdentityOperator(dim)

    concrete_ops = {
        op_name: op_lin,
        proj_name: proj_op,
        "1": id_op,
    }
    num_states = int(states.number_of_states)
    written: list[str] = []
    pair_arrays: list[tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]] = []

    for cur_op in (op_name, "1"):
        w_str = f"Tr[P:{proj_name} U(t)^+ {cur_op} U(t) delta(E-H) P:{proj_name}]"
        compiled = compile_word(
            w_str,
            basis={"name": "real", "operators": concrete_ops},
            params={
                "expansion": {
                    "moments": int(num_moments),
                    "states": num_states,
                    "window": window,
                },
                "time": {
                    "steps": int(ntime),
                    "tmax": float(tmax),
                    "bessel": str(bessel),
                },
                "reconstruction": {
                    "broadening_meV": float(broadening_meV),
                    "fermi_energies": [float(ef)],
                    "legacy_text_round_trip": bool(legacy_text_round_trip),
                },
            },
        )
        mu = run_plan(
            compiled.plan,
            ham,
            window,
            states,
            operators=compiled.operators,
            batch_size=int(num_moments),
            symmetrise=False,
            dt=compiled.dt,
        )
        mu.state_label = state_label
        stem = f"TimeProj{cur_op}-{proj_name}{label}KPM_M{int(num_moments)}_state{state_label}"
        mom_path = os.path.join(stage_dir, f"{stem}.chebmomTD")
        _core.WriteLegacyChebmomTD(mom_path, mu)
        mu_disk = _core.ReadChebmom(mom_path) if legacy_text_round_trip else mu

        times, vals = time_correlation_curve(
            mu_disk, float(broadening_meV), float(ef)
        )
        dat_name = f"mean{stem}EF{float(ef):.6f}JACKSON.dat"
        dat_path = os.path.join(
            out_dir if keep_top_level else stage_dir, dat_name
        )
        _core.WriteCurveText(dat_path, times, vals)
        if keep_top_level:
            written.append(dat_path)

        cap_path = os.path.join(captured_dir, f"{cap_idx:04d}_{dat_name}")
        shutil.copyfile(dat_path, cap_path)
        written.append(cap_path)
        cap_idx += 1

        loaded_t, loaded_v = np.loadtxt(dat_path, unpack=True)
        pair_arrays.append((loaded_t, loaded_v))

    return (
        pair_arrays[0][0],
        pair_arrays[0][1],
        pair_arrays[1][1],
        written,
        cap_idx,
    )


def run_tutorial06(
    out_dir: str,
    *,
    ey_l: int = 20,
    ey_ds: Sequence[float] = (0.05, 0.07, 0.10),
    ey_moments: int = 128,
    ey_steps: int = 80,
    ey_tmax: float = 260.0,
    dp_l: int = 20,
    dp_lambda: float = 0.2,
    dp_ws: Sequence[float] = (0.0, 6.0),
    dp_moments: int = 160,
    dp_steps: int = 120,
    dp_tmax: float = 45.0,
    ratelaw_l: int = 50,
    ratelaw_ds: Sequence[float] = (0.04, 0.06, 0.08, 0.10, 0.14),
    ratelaw_nr: int = 20,
    ratelaw_moments: int = 128,
    ratelaw_steps: int = 80,
    ratelaw_tmax: float = 260.0,
    crosscheck_l: int = 8,
    crosscheck_d: float = 0.10,
    crosscheck_nr: int = 20,
    crosscheck_moments: int = 96,
    crosscheck_steps: int = 60,
    crosscheck_tmax: float = 220.0,
    disorder_seed: int = 1,
    broadening_meV: float = 140.0,
    fermi_energies: Sequence[float] = (2.0,),
    bessel: str = "legacy",
    legacy_text_round_trip: bool = True,
    words: Sequence[str] = (
        "Tr[P:PZ U(t)^+ SZ U(t) delta(E-H) P:PZ]",
        "Tr[P:PZ U(t)^+ 1 U(t) delta(E-H) P:PZ]",
        "Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]",
        "Tr[P:PX U(t)^+ 1 U(t) delta(E-H) P:PX]",
    ),
) -> list[str]:
    """Run Tutorial 06 and write all top-level and captured ``.dat`` files into ``out_dir``."""
    for w_str in words:
        compile_word(w_str, params={"steps": int(ey_steps), "tmax": float(ey_tmax)})

    mod_ey, mod_rashba = _load_tutorial06_generators()
    os.makedirs(out_dir, exist_ok=True)
    captured_dir = os.path.join(out_dir, "captured")
    os.makedirs(captured_dir, exist_ok=True)

    ef = float(fermi_energies[0]) if fermi_energies else 2.0
    written: list[str] = []
    cap_idx = 0

    # 1. Part A: Elliott-Yafet exact trace on L = 20, D in (0.05, 0.07, 0.10) (captures 0000..0005)
    for d_val in ey_ds:
        with tempfile.TemporaryDirectory() as stage_dir:
            label = str(
                mod_ey.build_ey(
                    int(ey_l), float(d_val), int(disorder_seed), outdir=stage_dir
                )
            )
            dim = 2 * int(ey_l) * int(ey_l)
            states = _core.LocalStates(dim, list(range(dim)))
            _, _, _, w_paths, cap_idx = _eval_and_reconstruct_pair(
                stage_dir,
                label,
                "SZ",
                "PZ",
                int(ey_moments),
                int(ey_steps),
                float(ey_tmax),
                states,
                "exact",
                float(broadening_meV),
                ef,
                bessel,
                legacy_text_round_trip,
                out_dir,
                captured_dir,
                cap_idx,
            )
            written.extend(w_paths)

    # 2. Part B: Dyakonov-Perel exact trace on L = 20, lambda = 0.2, W in (0.0, 6.0) (captures 0006..0011)
    for w_val in dp_ws:
        with tempfile.TemporaryDirectory() as stage_dir:
            label = str(
                mod_rashba.build_rashba(
                    int(dp_l),
                    float(dp_lambda),
                    float(w_val),
                    int(disorder_seed),
                    outdir=stage_dir,
                )
            )
            dim = 2 * int(dp_l) * int(dp_l)
            states = _core.LocalStates(dim, list(range(dim)))
            _, _, _, w_paths, cap_idx = _eval_and_reconstruct_pair(
                stage_dir,
                label,
                "SZ",
                "PZ",
                int(dp_moments),
                int(dp_steps),
                float(dp_tmax),
                states,
                "exact",
                float(broadening_meV),
                ef,
                bessel,
                legacy_text_round_trip,
                out_dir,
                captured_dir,
                cap_idx,
            )
            written.extend(w_paths)
            if float(w_val) == float(dp_ws[-1]):
                _, _, _, w_paths_x, cap_idx = _eval_and_reconstruct_pair(
                    stage_dir,
                    label,
                    "SX",
                    "PX",
                    int(dp_moments),
                    int(dp_steps),
                    float(dp_tmax),
                    states,
                    "exact",
                    float(broadening_meV),
                    ef,
                    bessel,
                    legacy_text_round_trip,
                    out_dir,
                    captured_dir,
                    cap_idx,
                )
                written.extend(w_paths_x)

    # 3. rp_validate (1): EY rate-law random-phase Nr = 20 on L = 50 (captures 0012..0211)
    ratelaw_d_list: list[float] = []
    ratelaw_rates: list[float] = []
    for d_val in ratelaw_ds:
        with tempfile.TemporaryDirectory() as stage_dir:
            label = str(
                mod_ey.build_ey(
                    int(ratelaw_l),
                    float(d_val),
                    int(disorder_seed),
                    outdir=stage_dir,
                )
            )
            dim = 2 * int(ratelaw_l) * int(ratelaw_l)
            num_acc: npt.NDArray[np.float64] | None = None
            den_acc: npt.NDArray[np.float64] | None = None
            t_arr: npt.NDArray[np.float64] | None = None
            for r_idx in range(1, int(ratelaw_nr) + 1):
                states = _core.RandomPhaseStates(dim, 1, r_idx)
                tt, num_v, den_v, w_paths, cap_idx = _eval_and_reconstruct_pair(
                    stage_dir,
                    label,
                    "SZ",
                    "PZ",
                    int(ratelaw_moments),
                    int(ratelaw_steps),
                    float(ratelaw_tmax),
                    states,
                    states.label,
                    float(broadening_meV),
                    ef,
                    bessel,
                    legacy_text_round_trip,
                    out_dir,
                    captured_dir,
                    cap_idx,
                )
                written.extend(w_paths)
                num_acc = num_v if num_acc is None else num_acc + num_v
                den_acc = den_v if den_acc is None else den_acc + den_v
                t_arr = tt
            assert t_arr is not None and num_acc is not None and den_acc is not None
            good = t_arr > 0.0
            s_curve = num_acc[good] / float(np.mean(den_acc[good]))
            s_norm = s_curve / float(s_curve[0])
            t1_val = _fit_rate(t_arr[good], s_norm)
            ratelaw_d_list.append(float(d_val))
            ratelaw_rates.append(1.0 / t1_val)

    # 4. rp_validate (2): exact vs random-phase cross-check at L = 8, D = 0.10 (captures 0212..0253)
    with tempfile.TemporaryDirectory() as stage_dir:
        label_sm = str(
            mod_ey.build_ey(
                int(crosscheck_l),
                float(crosscheck_d),
                int(disorder_seed),
                outdir=stage_dir,
            )
        )
        dim_sm = 2 * int(crosscheck_l) * int(crosscheck_l)
        states_ex = _core.LocalStates(dim_sm, list(range(dim_sm)))
        _, _, _, w_paths_ex, cap_idx = _eval_and_reconstruct_pair(
            stage_dir,
            label_sm,
            "SZ",
            "PZ",
            int(crosscheck_moments),
            int(crosscheck_steps),
            float(crosscheck_tmax),
            states_ex,
            f"exact_{dim_sm}",
            float(broadening_meV),
            ef,
            bessel,
            legacy_text_round_trip,
            out_dir,
            captured_dir,
            cap_idx,
            keep_top_level=False,
        )
        written.extend(w_paths_ex)
        for r_idx in range(1, int(crosscheck_nr) + 1):
            states_rp = _core.RandomPhaseStates(dim_sm, 1, r_idx)
            _, _, _, w_paths_rp, cap_idx = _eval_and_reconstruct_pair(
                stage_dir,
                label_sm,
                "SZ",
                "PZ",
                int(crosscheck_moments),
                int(crosscheck_steps),
                float(crosscheck_tmax),
                states_rp,
                states_rp.label,
                float(broadening_meV),
                ef,
                bessel,
                legacy_text_round_trip,
                out_dir,
                captured_dir,
                cap_idx,
            )
            written.extend(w_paths_rp)

    # 5. Write ey_ratelaw.dat (and capture 0254_ey_ratelaw.dat)
    d_arr = np.array(ratelaw_d_list, dtype=np.float64)
    r_arr = np.array(ratelaw_rates, dtype=np.float64)
    p_fit, _ = np.polyfit(np.log(d_arr), np.log(r_arr), 1)
    ratelaw_path = os.path.join(out_dir, "ey_ratelaw.dat")
    np.savetxt(
        ratelaw_path,
        np.column_stack([d_arr, r_arr]),
        header=f"Dsf  1/T1  (p={float(p_fit):.3f})",
    )
    written.append(ratelaw_path)
    ratelaw_cap = os.path.join(captured_dir, f"{cap_idx:04d}_ey_ratelaw.dat")
    shutil.copyfile(ratelaw_path, ratelaw_cap)
    written.append(ratelaw_cap)

    return written


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``lsquant.apps.tutorial06``."""
    try:
        parser = argparse.ArgumentParser(prog="lsquant.apps.tutorial06")
        parser.add_argument("pos_out", nargs="?", default=None, help="Output directory")
        parser.add_argument(
            "--out", "-o", dest="out", default=None, help="Output directory"
        )
        args = parser.parse_args(argv)
        out_dir = (
            args.out
            if args.out is not None
            else (args.pos_out if args.pos_out is not None else ".")
        )
        paths = run_tutorial06(out_dir)
        for p in paths:
            print(p)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"tutorial06 error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
