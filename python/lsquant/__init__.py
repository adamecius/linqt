"""LSQUANT Python library for Kernel Polynomial Method calculations.

Exposes the C++ execution engine via :mod:`lsquant._core` alongside SymPy
spectral functions, damping kernels, moment plan builders, and rank-1/rank-2
curve reconstructions.
"""

from __future__ import annotations

import sys as _sys

try:
    from . import _core
except ImportError:
    import _core as _core  # type: ignore[no-redef]

    _sys.modules[f"{__name__}._core"] = _core

from . import kernels, moments, reconstruct, spectral_functions
from .kernels import jackson_damping, jackson_damping_2d
from .moments import dos_plan, kubo_plan, run_plan, window_from_gershgorin
from .reconstruct import (
    dos_curve,
    kubo_bastin,
    kubo_bastin_sea,
    kubo_bastin_surf,
    kubo_fft_grid,
    kubo_greenwood,
)
from .spectral_functions import (
    DeltaFunction,
    GreenRetarded,
    GreenRetardedDerivative,
    SpectralFunction,
    legacy_delta_table,
    legacy_green_tables,
)

__all__ = [
    "_core",
    "kernels",
    "moments",
    "reconstruct",
    "spectral_functions",
    "jackson_damping",
    "jackson_damping_2d",
    "dos_plan",
    "kubo_plan",
    "window_from_gershgorin",
    "run_plan",
    "dos_curve",
    "kubo_greenwood",
    "kubo_bastin",
    "kubo_bastin_sea",
    "kubo_bastin_surf",
    "kubo_fft_grid",
    "SpectralFunction",
    "DeltaFunction",
    "GreenRetarded",
    "GreenRetardedDerivative",
    "legacy_delta_table",
    "legacy_green_tables",
]
