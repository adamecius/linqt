"""LSQUANT Python library for Kernel Polynomial Method calculations.

Exposes the C++ execution engine via :mod:`lsquant._core` alongside SymPy
spectral functions, damping kernels, moment plan builders, and rank-1 curve
reconstruction.
"""

from __future__ import annotations

import sys as _sys

try:
    from . import _core
except ImportError:
    import _core as _core  # type: ignore[no-redef]

    _sys.modules[f"{__name__}._core"] = _core

from . import kernels, moments, reconstruct, spectral_functions
from .kernels import jackson_damping
from .moments import dos_plan, kubo_plan, run_plan
from .reconstruct import dos_curve
from .spectral_functions import DeltaFunction, SpectralFunction, legacy_delta_table

__all__ = [
    "_core",
    "kernels",
    "moments",
    "reconstruct",
    "spectral_functions",
    "jackson_damping",
    "dos_plan",
    "kubo_plan",
    "run_plan",
    "dos_curve",
    "SpectralFunction",
    "DeltaFunction",
    "legacy_delta_table",
]
