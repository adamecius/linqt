"""Chebyshev truncation damping kernels.

Provides Jackson kernel damping factors computed with standard library ``math``
trigonometric functions in the exact legacy operation order to guarantee
bit-identical moment damping.
"""

from __future__ import annotations

import math
import numpy as np
import numpy.typing as npt


def jackson_damping(num_moments: int) -> npt.NDArray[np.float64]:
    """Compute the Jackson kernel damping weights ``g_m`` for ``m = 0 .. M - 1``.

    Solves the problem of suppressing Gibbs oscillations in a truncated
    Chebyshev series of length ``M = num_moments`` while preserving non-negativity
    of spectral densities. Evaluates the closed-form Jackson weights using
    Python's ``math`` module in the exact operation order of legacy
    ``chebyshev::Moments::JacksonKernel``:

    .. math::

        \\phi_J = \\frac{\\pi}{M + 1}, \\qquad
        g_m = \\left((M - m + 1)\\cos(m\\phi_J)
              + \\frac{\\sin(m\\phi_J)}{\\tan(\\phi_J)}\\right)\\frac{\\phi_J}{\\pi}.

    All quantities are dimensionless and lie in ``(0, 1]`` with ``g_0 = 1``.
    Time and space complexity are ``O(num_moments)``.

    Parameters
    ----------
    num_moments : int
        Effective number of Chebyshev moments ``M > 0``.

    Returns
    -------
    numpy.ndarray
        1D ``float64`` array of length ``num_moments`` containing ``g_m``.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments <= 0``.

    Notes
    -----
    Pitfall: Do not replace ``math.cos``, ``math.sin``, and ``math.tan`` with
    vectorised NumPy ufuncs; SIMD trigonometric implementations can differ by
    1 ULP and break byte-for-byte golden curve reproduction.
    """
    if num_moments <= 0:
        raise ValueError("jackson_damping: num_moments must be positive")

    mom = float(num_moments)
    phi = math.pi / (mom + 1.0)
    weights = np.empty(num_moments, dtype=np.float64)
    for m_idx in range(num_moments):
        m = float(m_idx)
        weights[m_idx] = (
            ((mom - m + 1.0) * math.cos(phi * m) + math.sin(phi * m) / math.tan(phi))
            * phi
            / math.pi
        )
    return weights
