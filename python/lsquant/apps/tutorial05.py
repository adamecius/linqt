"""Reproduce Tutorial 05 (diffusive charge transport on a 2D square lattice) through the Python API.

Imports ``examples/05_charge_transport_square/make_square.py``, generates the
``L = 200`` (``N = 40000``) disordered square-lattice operators for
``W in (0.5, 1.0, 2.0, 4.0)`` eV with disorder seed ``1`` and ``KPM_SEED = 1``,
computes the MSD moments (``M = 256``, ``160`` steps, ``tmax = 150`` fs) and
reconstructs at ``E_F = 2.0`` eV (``80`` meV broadening), and computes the DOS
moments (``M = 1024``) and reconstructs at ``80`` meV broadening.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
import tempfile
from types import ModuleType

from .. import _core
from ..moments import (
    dos_plan,
    msd_plan,
    run_plan,
    time_evolution_series,
    time_step_from_tmax,
)
from ..reconstruct import dos_curve, time_correlation_curve

L_SIDE: int = 200
W_LIST: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0)
DISORDER_SEED: int = 1
KPM_SEED: int = 1
M_MSD: int = 256
NTIME: int = 160
TMAX: float = 150.0
EF_TARGET: float = 2.0
M_DOS: int = 1024
BROAD_MEV: float = 80.0


def _load_make_square() -> ModuleType:
    """Import ``examples/05_charge_transport_square/make_square.py``."""
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..")
    )
    script_path = os.path.join(
        repo_root, "examples", "05_charge_transport_square", "make_square.py"
    )
    if not os.path.isfile(script_path):
        raise FileNotFoundError(
            f"tutorial05: make_square.py not found at {script_path}"
        )
    spec = importlib.util.spec_from_file_location(
        "lsquant_tutorial05_make_square", script_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"tutorial05: could not load spec for {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_bounds(bounds_path: str) -> _core.SpectralWindow:
    """Read ``lo hi`` from a legacy ``BOUNDS`` file."""
    with open(bounds_path, encoding="utf-8") as fh:
        tokens = fh.read().split()
    return _core.SpectralWindow(float(tokens[0]), float(tokens[1]))


def run_tutorial05(out_dir: str) -> list[str]:
    """Run Tutorial 05 and write the eight MSD and DOS ``.dat`` curves into ``out_dir``."""
    make_square = _load_make_square()
    os.makedirs(out_dir, exist_ok=True)
    captured_dir = os.path.join(out_dir, "captured")
    os.makedirs(captured_dir, exist_ok=True)

    written: list[str] = []
    cap_idx = 0
    dt = time_step_from_tmax(NTIME, TMAX)
    tmax_str = f"{TMAX:g}"

    for w_val in W_LIST:
        with tempfile.TemporaryDirectory() as stage_dir:
            label = str(
                make_square.build_square(L_SIDE, w_val, DISORDER_SEED, outdir=stage_dir)
            )
            op_dir = os.path.join(stage_dir, "operators")
            window = _read_bounds(os.path.join(stage_dir, "BOUNDS"))

            ham = _core.ReadCsrText(os.path.join(op_dir, f"{label}.HAM.CSR"))
            vx_mat = _core.ReadCsrText(os.path.join(op_dir, f"{label}.VX.CSR"))
            vx_op = _core.RescaledOperator(vx_mat, window)
            dim = ham.dimension

            # 1. MSD at E_F = 2.0 eV
            states_msd = _core.RandomPhaseStates(dim, 1, KPM_SEED)
            # Byte identity with the legacy captures requires the GSL Bessel values.
            coeffs, phase = time_evolution_series(window, dt, use_legacy_bessel=True)
            plan_msd = msd_plan(
                M_MSD,
                NTIME,
                1,
                coeffs,
                velocity_index=0,
                phase=phase,
                word="[VX, U(t)]^+ delta [VX, U(t)]",
            )
            mu_msd = run_plan(
                plan_msd,
                ham,
                window,
                states_msd,
                operators=[vx_op],
                batch_size=M_MSD,
                symmetrise=False,
                dt=dt,
            )
            msd_stem = (
                f"CorrelationVX{label}KPM_M{M_MSD}_state{tmax_str}{states_msd.label}"
            )
            msd_mom_path = os.path.join(stage_dir, f"{msd_stem}.chebmomTD")
            _core.WriteLegacyChebmomTD(msd_mom_path, mu_msd)
            mu_msd_disk = _core.ReadChebmom(msd_mom_path)

            times, msd_vals = time_correlation_curve(
                mu_msd_disk, BROAD_MEV, EF_TARGET
            )
            msd_dat_name = f"mean{msd_stem}EF{EF_TARGET:.6f}JACKSON.dat"
            msd_dat_path = os.path.join(out_dir, msd_dat_name)
            _core.WriteCurveText(msd_dat_path, times, msd_vals)
            written.append(msd_dat_path)

            msd_cap_path = os.path.join(
                captured_dir, f"{cap_idx:04d}_{msd_dat_name}"
            )
            shutil.copyfile(msd_dat_path, msd_cap_path)
            cap_idx += 1

            # 2. DOS at M = 1024
            states_dos = _core.RandomPhaseStates(dim, 1, KPM_SEED)
            plan_dos = dos_plan(M_DOS, 1, word="1 delta")
            mu_dos = run_plan(
                plan_dos, ham, window, states_dos, batch_size=M_DOS
            )
            dos_stem = f"SpectralOp1{label}KPM_M{M_DOS}_state{states_dos.label}"
            dos_mom_path = os.path.join(stage_dir, f"{dos_stem}.chebmom1D")
            _core.WriteLegacyChebmom1D(dos_mom_path, mu_dos)
            mu_dos_disk = _core.ReadChebmom(dos_mom_path)

            energies, dos_vals = dos_curve(mu_dos_disk, BROAD_MEV)
            dos_dat_name = f"mean{dos_stem}JACKSON.dat"
            dos_dat_path = os.path.join(out_dir, dos_dat_name)
            _core.WriteCurveText(dos_dat_path, energies, dos_vals)
            written.append(dos_dat_path)

            dos_cap_path = os.path.join(
                captured_dir, f"{cap_idx:04d}_{dos_dat_name}"
            )
            shutil.copyfile(dos_dat_path, dos_cap_path)
            cap_idx += 1

    return written


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``lsquant.apps.tutorial05``."""
    try:
        parser = argparse.ArgumentParser(prog="lsquant.apps.tutorial05")
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
        paths = run_tutorial05(out_dir)
        for p in paths:
            print(p)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"tutorial05 error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
