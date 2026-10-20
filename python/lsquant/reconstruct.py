"""Rank-1 Chebyshev spectral reconstruction routines.

Reconstructs the extensive density of states curve from a rank-1
:class:`~lsquant._core.ChebyshevMoments` tensor using Jackson damping and the
legacy uniform reconstruction grid.
"""

from __future__ import annotations

import math
import numpy as np
import numpy.typing as npt

from . import _core
from .kernels import jackson_damping
from .spectral_functions import legacy_delta_table


def dos_curve(
    moments: _core.ChebyshevMoments,
    broadening_meV: float,
    alpha: float = 0.95,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Reconstruct the extensive density of states curve from rank-1 moments.

    Solves the problem of turning a rank-1 Chebyshev moment vector ``mu_m`` and
    a target energy resolution ``broadening_meV`` into the physical energy grid
    ``E`` (in eV) and extensive density of states ``N * rho(E)`` (in 1/eV),
    reproducing the operation order of legacy ``spectralFunctionFromChebmom``
    byte for byte:

    1. ``eta = 2.0 * broadening_meV / 1000 / BandWidth``,
       ``M_eff = min(ceil(pi / eta), M)``.
    2. Truncate moments to ``M_eff``, multiply complex ``mu[m]`` by the Jackson
       kernel ``g(m, M_eff)``, then take the real part ``mu_real[m] = mu[m].real``.
    3. Build ``num_div = 30 * M_eff`` uniform points
       ``x_i = -alpha + i * (2 * alpha) / (num_div - 1)``.
    4. Contract ``rho_i = sum_{m=0}^{M_eff-1} delta_chebF(x_i, m) * mu_real[m]``
       sequentially in ``m`` via :func:`lsquant._core.ContractRank1`.
    5. Map to physical energy ``E_i = x_i / Scale + Shift`` and extensive value
       ``value_i = rho_i * (SystemSize * Scale)``.

    Complexity is ``O(M_eff * num_div) = O(30 * M_eff^2)`` time and memory.

    Parameters
    ----------
    moments : lsquant._core.ChebyshevMoments
        Rank-1 Chebyshev moment tensor with populated ``window`` and
        ``system_size``.
    broadening_meV : float
        Target Jackson Gaussian broadening in meV (``broadening_meV > 0``).
    alpha : float, optional
        Dimensionless Chebyshev domain cutoff keeping the grid inside
        ``[-alpha, alpha]`` away from ``x = +-1`` (default ``0.95``).

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Pair ``(energies, values)`` of 1D ``float64`` arrays of length
        ``30 * M_eff``, where ``energies`` is in eV and ``values`` is the
        extensive density of states ``SystemSize * rho(E)`` in 1/eV.

    Warnings
    --------
    Raises ``ValueError`` if ``moments.rank != 1``, ``broadening_meV <= 0``,
    ``!(0 < alpha < 1)``, or the spectral window has non-positive width.

    Notes
    -----
    Pitfall: The returned ``values`` array is the extensive density of states
    ``N * rho(E)`` written by ``spectralFunctionFromChebmom``; divide by
    ``moments.system_size`` when plotting the intensive per-site density of
    states.
    """
    if moments.rank != 1:
        raise ValueError("dos_curve: moments must have rank 1")
    if broadening_meV <= 0.0:
        raise ValueError("dos_curve: broadening_meV must be positive")
    if not (0.0 < alpha < 1.0):
        raise ValueError("dos_curve: alpha must lie strictly between 0 and 1")

    band_width = (moments.window.upper - moments.window.lower) * 1.0
    if band_width <= 0.0:
        raise ValueError("dos_curve: spectral window must have positive width")

    eta = 2.0 * float(broadening_meV) / 1000.0 / band_width
    max_mom = int(math.ceil(math.pi / eta))
    total_mom = int(moments.dims[0])
    if max_mom > total_mom:
        max_mom = total_mom

    g = jackson_damping(max_mom)
    damped = np.empty(max_mom, dtype=np.complex128)
    raw = moments.data
    for m_idx in range(max_mom):
        damped[m_idx] = raw[m_idx] * g[m_idx]

    num_div = 30 * max_mom
    x_grid = [
        -alpha + i * (2.0 * alpha) / (num_div - 1) for i in range(num_div)
    ]

    table = legacy_delta_table(x_grid, max_mom)
    rho = _core.ContractRank1(damped.real, table)

    scale = moments.window.scale
    shift = moments.window.shift
    prefactor = float(moments.system_size) * scale

    energies = np.empty(num_div, dtype=np.float64)
    values = np.empty(num_div, dtype=np.float64)
    for i in range(num_div):
        energies[i] = x_grid[i] / scale + shift
        values[i] = float(rho[i]) * prefactor

    return energies, values
