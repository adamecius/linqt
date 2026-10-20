"""Reproduce Tutorial 01 (1D chain density of states) through the Python API.

Builds the periodic 1D tight-binding chains for ``N = 256, 1024, 4096`` by
importing ``examples/01_dos_1d_chain/make_chain.py``, evaluates ``M = 512``
exact-trace Chebyshev moments through :class:`~lsquant._core.MomentEngine`,
writes the legacy ``.chebmom1D`` files, and reconstructs the density of states
at ``20, 80, 24, 10`` meV (as well as the README ``120`` meV step) with
:func:`~lsquant.reconstruct.dos_curve`.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
from types import ModuleType

from .. import _core
from ..moments import dos_plan, run_plan
from ..reconstruct import dos_curve


def _load_make_chain() -> ModuleType:
    """Import ``examples/01_dos_1d_chain/make_chain.py`` from the repository tree.

    Returns
    -------
    types.ModuleType
        Loaded ``make_chain`` module.

    Warnings
    --------
    Raises ``FileNotFoundError`` if ``examples/01_dos_1d_chain/make_chain.py``
    cannot be located relative to the package path.

    Notes
    -----
    Pitfall: Loads ``make_chain.py`` via :mod:`importlib.util` so that no copy
    of the model generator exists in the library.
    """
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..")
    )
    script_path = os.path.join(repo_root, "examples", "01_dos_1d_chain", "make_chain.py")
    if not os.path.isfile(script_path):
        raise FileNotFoundError(f"tutorial01: make_chain.py not found at {script_path}")

    spec = importlib.util.spec_from_file_location("lsquant_tutorial01_make_chain", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"tutorial01: could not load spec for {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_tutorial01(out_dir: str) -> list[str]:
    """Generate the Tutorial 01 operators, exact moments, and DOS curves in ``out_dir``.

    Solves the problem of reproducing every ``.chebmom1D`` and ``JACKSON.dat``
    file captured from ``examples/01_dos_1d_chain`` using the new Python API.
    Writes into ``out_dir`` both the unprefixed files and the numbered sequence
    matching ``goldens/stage0/tutorials/01/captured/`` (in ``out_dir`` and
    ``out_dir/captured/``).

    Parameters
    ----------
    out_dir : str
        Destination directory for generated operators, moment files, and curves.

    Returns
    -------
    list[str]
        Ordered list of primary output file paths written.

    Warnings
    --------
    Overwrites existing tutorial output files in ``out_dir``.

    Notes
    -----
    Pitfall: Because ``lsqplot.py`` reconstructs ``moments_chain1d_N256_M512.chebmom1D``
    at four broadenings (``20, 80, 24, 10`` meV) in separate temporary directories,
    all four reconstructions share the legacy basename
    ``meanmoments_chain1d_N256_M512JACKSON.dat`` and are disambiguated by their
    capture prefixes ``0013_``, ``0019_``, ``0020_``, and ``0021_`` (as well as
    ``eta<broadening>/`` subdirectories).
    """
    make_chain = _load_make_chain()
    os.makedirs(out_dir, exist_ok=True)
    op_dir = os.path.join(out_dir, "operators")
    cap_dir = os.path.join(out_dir, "captured")
    os.makedirs(op_dir, exist_ok=True)
    os.makedirs(cap_dir, exist_ok=True)

    sizes = (256, 1024, 4096)
    num_moments = 512
    hopping_t = -1.0
    window = _core.SpectralWindow(-2.0, 2.0)

    moments_by_size: dict[int, _core.ChebyshevMoments] = {}
    spectral_mom_paths: dict[int, str] = {}
    short_mom_paths: dict[int, str] = {}
    written: list[str] = []

    for n_sites in sizes:
        label = f"chain1d_N{n_sites}"
        csr_path = os.path.join(op_dir, f"{label}.HAM.CSR")
        desc_path = os.path.join(op_dir, f"{label}.HAM.desc")
        bounds_path = os.path.join(out_dir, "BOUNDS")
        state_path = os.path.join(out_dir, f"exact_{label}")

        make_chain.write_csr(csr_path, n_sites, hopping_t)
        make_chain.write_desc(desc_path, label, n_sites, hopping_t, csr_path=csr_path)
        make_chain.write_bounds(bounds_path, hopping_t)
        make_chain.write_state(state_path, n_sites)

        ham = _core.ReadCsrText(csr_path)
        states = _core.LocalStates(n_sites, list(range(n_sites)))
        plan = dos_plan(num_moments, n_sites)
        mu = run_plan(plan, ham, window, states, batch_size=num_moments)
        mu.state_label = f"exact_{label}"
        moments_by_size[n_sites] = mu

        spec_name = f"SpectralOp1{label}KPM_M{num_moments}_stateexact_{label}.chebmom1D"
        short_name = f"moments_{label}_M{num_moments}.chebmom1D"
        spec_path = os.path.join(out_dir, spec_name)
        short_path = os.path.join(out_dir, short_name)

        _core.WriteLegacyChebmom1D(spec_path, mu)
        shutil.copyfile(spec_path, short_path)

        spectral_mom_paths[n_sites] = spec_path
        short_mom_paths[n_sites] = short_path
        written.append(spec_path)
        written.append(short_path)

    # README step 4: 120 meV reconstruction on N=256
    e_120, v_120 = dos_curve(moments_by_size[256], 120.0)
    readme_dat = os.path.join(
        out_dir,
        "meanSpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256JACKSON.dat",
    )
    _core.WriteCurveText(readme_dat, e_120, v_120)
    written.append(readme_dat)

    # Capture sequence 0000..0011: figure_moments and figure_dos moment reads
    capture_moment_entries = [
        (0, 256, True),
        (1, 256, False),
        (2, 1024, True),
        (3, 1024, False),
        (4, 4096, True),
        (5, 4096, False),
        (6, 256, True),
        (7, 256, False),
        (8, 1024, True),
        (9, 1024, False),
        (10, 4096, True),
        (11, 4096, False),
        (12, 256, True),
        (14, 1024, True),
        (16, 4096, True),
        (18, 256, True),
    ]
    for idx, n_sites, is_spec in capture_moment_entries:
        src = spectral_mom_paths[n_sites] if is_spec else short_mom_paths[n_sites]
        cname = f"{idx:04d}_{os.path.basename(src)}"
        shutil.copyfile(src, os.path.join(cap_dir, cname))

    # Panel (a): fixed broadening 20 meV across N = 256, 1024, 4096
    panel_a_Specs = [
        (13, 256, 20.0),
        (15, 1024, 20.0),
        (17, 4096, 20.0),
    ]
    for idx, n_sites, eta in panel_a_Specs:
        energies, values = dos_curve(moments_by_size[n_sites], eta)
        base_name = f"meanmoments_chain1d_N{n_sites}_M{num_moments}JACKSON.dat"
        prefixed_name = f"{idx:04d}_{base_name}"
        prefixed_path = os.path.join(out_dir, prefixed_name)
        _core.WriteCurveText(prefixed_path, energies, values)
        shutil.copyfile(prefixed_path, os.path.join(cap_dir, prefixed_name))

        eta_dir = os.path.join(out_dir, f"eta{int(eta)}")
        os.makedirs(eta_dir, exist_ok=True)
        shutil.copyfile(prefixed_path, os.path.join(eta_dir, base_name))
        if n_sites in (1024, 4096):
            unprefixed_path = os.path.join(out_dir, base_name)
            shutil.copyfile(prefixed_path, unprefixed_path)
        written.append(prefixed_path)

    # Panel (b): N = 256 at broadenings 80, 24, 10 meV
    panel_b_Specs = [
        (19, 256, 80.0),
        (20, 256, 24.0),
        (21, 256, 10.0),
    ]
    for idx, n_sites, eta in panel_b_Specs:
        energies, values = dos_curve(moments_by_size[n_sites], eta)
        base_name = f"meanmoments_chain1d_N{n_sites}_M{num_moments}JACKSON.dat"
        prefixed_name = f"{idx:04d}_{base_name}"
        prefixed_path = os.path.join(out_dir, prefixed_name)
        _core.WriteCurveText(prefixed_path, energies, values)
        shutil.copyfile(prefixed_path, os.path.join(cap_dir, prefixed_name))

        eta_dir = os.path.join(out_dir, f"eta{int(eta)}")
        os.makedirs(eta_dir, exist_ok=True)
        shutil.copyfile(prefixed_path, os.path.join(eta_dir, base_name))
        written.append(prefixed_path)

    return written


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``lsquant.apps.tutorial01``.

    Parameters
    ----------
    argv : list[str] or None, optional
        Command-line arguments.

    Returns
    -------
    int
        Exit status (``0`` on success, ``1`` on error).

    Warnings
    --------
    Prints any uncaught exception to ``stderr`` and returns ``1``.

    Notes
    -----
    Pitfall: Accepts both ``--out <dir>`` and a positional ``<dir>`` argument,
    defaulting to the current working directory ``.``.
    """
    try:
        parser = argparse.ArgumentParser(prog="lsquant.apps.tutorial01")
        parser.add_argument("pos_out", nargs="?", default=None, help="Output directory")
        parser.add_argument("--out", "-o", dest="out", default=None, help="Output directory")
        args = parser.parse_args(argv)
        out_dir = args.out if args.out is not None else (args.pos_out if args.pos_out is not None else ".")
        paths = run_tutorial01(out_dir)
        for p in paths:
            print(p)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"tutorial01 error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
