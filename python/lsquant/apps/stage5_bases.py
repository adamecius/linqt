"""Stage-5 gate program evaluating DOS and Kubo words across Real, Bloch, and Eigen bases.

Invocation contract:

.. code-block:: bash

   python -m lsquant.apps.stage5_bases \\
     --model-dir <dir with <seed>_hr.dat, <seed>.uc, <seed>.xyz> \\
     --seed graphene \\
     --basis real|bloch|eigen \\
     --supercell 32 32 1 \\
     --disorder-w 0.0 \\
     --word dos|kuboxx \\
     --moments 256 \\
     --trace exact|random --states 10 --seed 12345 \\
     --window -10 10 \\
     --broadening-meV 100 \\
     --out <dir>

In all three bases, ``stage5_bases`` runs the identical :class:`~lsquant._core.Plan`
(:func:`~lsquant.moments.dos_plan` or :func:`~lsquant.moments.kubo_plan`) through
:class:`~lsquant._core.MomentEngine` and reconstructs with
:func:`~lsquant.reconstruct.dos_curve` or
:func:`~lsquant.reconstruct.kubo_greenwood`, writing ``.chebmom`` and ``.dat``
to ``--out`` and printing the written paths to ``stdout``.
"""

from __future__ import annotations

import argparse
import os
import sys

from .. import _core
from ..bases import make_basis
from ..models.wannier import load_wannier_model
from ..moments import dos_plan, kubo_plan, run_plan
from ..reconstruct import dos_curve, kubo_greenwood


def _resolve_seeds(
    model_dir: str,
    seed_tokens: list[str] | None,
    explicit_label: str | None,
    explicit_rng_seed: int | None,
) -> tuple[str | None, int]:
    """Disambiguate ``--seed <model_name>`` and ``--seed <int_rng_seed>`` when both are passed."""
    model_seed: str | None = explicit_label
    rng_seed: int | None = explicit_rng_seed

    for tok in seed_tokens or []:
        hr_candidate = os.path.join(model_dir, f"{tok}_hr.dat")
        if os.path.isfile(hr_candidate):
            model_seed = tok
            continue
        try:
            val = int(tok)
            rng_seed = val
        except ValueError:
            model_seed = tok

    if rng_seed is None:
        rng_seed = 12345
    return model_seed, int(rng_seed)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse and validate command-line options for ``stage5_bases``."""
    parser = argparse.ArgumentParser(prog="lsquant.apps.stage5_bases")
    parser.add_argument(
        "--model-dir",
        required=True,
        dest="model_dir",
        help="Directory containing <seed>_hr.dat, <seed>.uc, <seed>.xyz",
    )
    parser.add_argument(
        "--seed",
        action="append",
        default=[],
        help="Model seed prefix (e.g. 'graphene') and/or PRNG integer seed (e.g. 12345)",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Explicit model seed prefix if --seed is used only for the PRNG seed",
    )
    parser.add_argument(
        "--rng-seed",
        type=int,
        default=None,
        dest="rng_seed",
        help="Explicit integer PRNG seed",
    )
    parser.add_argument(
        "--basis",
        required=True,
        choices=("real", "bloch", "eigen"),
        help="Basis adapter ('real', 'bloch', or 'eigen')",
    )
    parser.add_argument(
        "--supercell",
        nargs=3,
        type=int,
        required=True,
        metavar=("N1", "N2", "N3"),
        help="Supercell dimensions N1 N2 N3",
    )
    parser.add_argument(
        "--disorder-w",
        type=float,
        default=0.0,
        dest="disorder_w",
        help="Real-space Anderson box disorder strength W in eV (default 0.0)",
    )
    parser.add_argument(
        "--word",
        required=True,
        choices=("dos", "kuboxx", "kuboxy"),
        help="Observable word ('dos', 'kuboxx', or 'kuboxy')",
    )
    parser.add_argument(
        "--moments",
        type=int,
        required=True,
        help="Number of Chebyshev moments M",
    )
    parser.add_argument(
        "--trace",
        choices=("exact", "random"),
        default="exact",
        help="Trace evaluation mode ('exact' or 'random')",
    )
    parser.add_argument(
        "--states",
        type=int,
        default=None,
        help="Number of random-phase trace vectors (when --trace random)",
    )
    parser.add_argument(
        "--window",
        nargs="+",
        required=True,
        help="Spectral window '<lower> <upper>' in eV or 'gershgorin'",
    )
    parser.add_argument(
        "--broadening-meV",
        type=float,
        required=True,
        dest="broadening_meV",
        help="Jackson Gaussian broadening in meV",
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=None,
        help="Chebyshev batch size (default M)",
    )
    parser.add_argument(
        "--out",
        default=".",
        help="Output directory",
    )

    args = parser.parse_args(argv)
    if args.moments <= 0:
        raise ValueError("stage5_bases: --moments must be positive")
    if args.broadening_meV <= 0.0:
        raise ValueError("stage5_bases: --broadening-meV must be positive")
    if args.disorder_w < 0.0:
        raise ValueError("stage5_bases: --disorder-w must be non-negative")
    return args


def _resolve_window(
    window_tokens: list[str], basis_obj: object
) -> _core.SpectralWindow:
    """Resolve ``--window <lo> <hi>`` or ``--window gershgorin``."""
    if len(window_tokens) == 1 and window_tokens[0].strip().lower() == "gershgorin":
        return basis_obj.window()  # type: ignore[attr-defined]
    if len(window_tokens) == 2:
        lo = float(window_tokens[0])
        hi = float(window_tokens[1])
        if not (hi > lo):
            raise ValueError("stage5_bases: --window upper bound must exceed lower bound")
        return _core.SpectralWindow(lo, hi)
    raise ValueError(
        "stage5_bases: --window must be either '<lo> <hi>' or 'gershgorin'"
    )


def main(argv: list[str] | None = None) -> int:
    """Run ``stage5_bases`` from command-line arguments."""
    try:
        args = _parse_args(argv)
        os.makedirs(args.out, exist_ok=True)

        model_seed, rng_seed = _resolve_seeds(
            args.model_dir, args.seed, args.label, args.rng_seed
        )
        model = load_wannier_model(args.model_dir, seed=model_seed)
        basis = make_basis(
            args.basis,
            model,
            args.supercell,
            disorder_w=float(args.disorder_w),
            seed=rng_seed,
        )

        window = _resolve_window(args.window, basis)
        state_gen = basis.states(
            trace_mode=args.trace,
            num_states=args.states,
            seed=rng_seed,
        )

        num_moments = int(args.moments)
        batch_size = int(args.batch) if args.batch is not None else num_moments
        word_key = args.word.strip().lower()

        if word_key == "dos":
            plan = dos_plan(
                num_moments,
                int(state_gen.number_of_states),
                word="1 delta",
            )
            mu = run_plan(
                plan,
                basis.hamiltonian(),
                window,
                state_gen,
                operators=[],
                batch_size=batch_size,
                symmetrise=False,
            )
            if args.trace == "exact":
                mu.number_of_states = basis.dimension
            energies, values = dos_curve(mu, float(args.broadening_meV))
            stem = (
                f"SpectralOp1{model.seed_name}_{args.basis}_KPM_M{num_moments}"
                f"_state{mu.state_label}"
            )
            cheb_path = os.path.join(args.out, f"{stem}.chebmom")
            dat_path = os.path.join(args.out, f"mean_{stem}JACKSON.dat")
        elif word_key in ("kuboxx", "kuboxy"):
            op_vx = basis.operator("VX")
            if word_key == "kuboxx":
                operators = [op_vx]
                left_idx, right_idx = 0, 0
                word_str = "VX delta VX delta"
                op_tag = "VX-VX"
            else:
                op_vy = basis.operator("VY")
                operators = [op_vx, op_vy]
                left_idx, right_idx = 0, 1
                word_str = "VX delta VY delta"
                op_tag = "VX-VY"
            plan = kubo_plan(
                num_moments,
                left_index=left_idx,
                right_index=right_idx,
                word=word_str,
            )
            mu = run_plan(
                plan,
                basis.hamiltonian(),
                window,
                state_gen,
                operators=operators,
                batch_size=batch_size,
                symmetrise=(word_key == "kuboxx"),
            )
            if args.trace == "exact":
                mu.number_of_states = basis.dimension
            energies, values = kubo_greenwood(mu, float(args.broadening_meV))
            stem = (
                f"NonEqOp{op_tag}{model.seed_name}_{args.basis}_KPM_M"
                f"{num_moments}x{num_moments}_state{mu.state_label}"
            )
            cheb_path = os.path.join(args.out, f"{stem}.chebmom")
            dat_path = os.path.join(args.out, f"KuboGreenWood_{stem}JACKSON.dat")
        else:
            raise ValueError(f"stage5_bases: unsupported word '{args.word}'")

        _core.WriteChebmom(cheb_path, mu)
        _core.WriteCurveText(dat_path, energies, values)
        print(cheb_path)
        print(dat_path)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"stage5_bases error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
