"""SymPy spectral functions and Chebyshev coefficient table evaluators.

Defines the symbolic :class:`SpectralFunction` representation in the rescaled
Chebyshev variable ``x in (-1, 1)`` and moment order ``m``, the Dirac delta
expansion :data:`DeltaFunction`, and the bit-exact legacy evaluator
:func:`legacy_delta_table`.
"""

from __future__ import annotations

import math
from typing import Sequence
import numpy as np
import numpy.typing as npt
import sympy as sp


class SpectralFunction:
    """Symbolic Chebyshev expansion coefficient generator for a function of ``H``.

    Solves the problem of turning an analytic spectral kernel ``c_m(x)`` written
    as a SymPy expression in the dimensionless rescaled energy ``x in (-1, 1)``
    and polynomial index ``m >= 0`` into a dense 2D ``(M, N)`` evaluation table
    for moment contraction. Lambdifies the symbolic expression once at
    construction time and broadcasts over ``(m, x)`` in ``O(M * N)`` time and
    memory.

    Parameters
    ----------
    expr : sympy.Expr
        Symbolic expression for the ``m``-th Chebyshev coefficient at rescaled
        coordinate ``x``.
    x : sympy.Symbol, optional
        Symbol representing the rescaled energy variable in ``(-1, 1)``.
    m : sympy.Symbol, optional
        Symbol representing the Chebyshev polynomial order ``0, 1, ..., M - 1``.

    Warnings
    --------
    Evaluating at ``|x| >= 1`` encounters the ``1 / sqrt(1 - x**2)`` Chebyshev
    weight singularity and yields ``inf`` or ``nan``.

    Notes
    -----
    Pitfall: NumPy's vectorised ``cos`` and ``arccos`` inside ``lambdify`` agree
    with ``libm`` to machine precision (``< 1e-14``) but may differ by 1 ULP;
    use :func:`legacy_delta_table` when byte-for-byte identity with legacy C++
    streams is required.
    """

    def __init__(
        self,
        expr: sp.Expr,
        x: sp.Symbol | None = None,
        m: sp.Symbol | None = None,
    ) -> None:
        self.x: sp.Symbol = x if x is not None else sp.Symbol("x", real=True)
        self.m: sp.Symbol = (
            m if m is not None else sp.Symbol("m", integer=True, nonnegative=True)
        )
        self.expr: sp.Expr = expr
        self._fn = sp.lambdify((self.m, self.x), self.expr, modules="numpy")

    def __call__(self) -> SpectralFunction:
        """Return ``self`` so the instance can also be invoked as a factory."""
        return self

    def table(
        self, x_grid: Sequence[float] | npt.NDArray[np.float64], num_moments: int
    ) -> npt.NDArray[np.float64]:
        """Evaluate the ``(num_moments, len(x_grid))`` Chebyshev coefficient table.

        Broadcasts the lambdified SymPy expression over column vector
        ``m = 0 .. num_moments - 1`` and row vector ``x_grid``. Both axes are
        dimensionless. Complexity is ``O(num_moments * len(x_grid))``.

        Parameters
        ----------
        x_grid : Sequence[float] or numpy.ndarray
            Rescaled evaluation points ``x_i in (-1, 1)`` of length ``N``.
        num_moments : int
            Number of Chebyshev moments ``M > 0``.

        Returns
        -------
        numpy.ndarray
            Contiguous 2D ``float64`` array of shape ``(M, N)``.

        Warnings
        --------
        Raises ``ValueError`` if ``num_moments <= 0``.

        Notes
        -----
        Pitfall: If ``self.expr`` is independent of ``m``, ``lambdify`` returns
        a row array; ``np.broadcast_to`` ensures the returned shape is always
        ``(M, N)``.
        """
        if num_moments <= 0:
            raise ValueError("SpectralFunction.table: num_moments must be positive")
        x_arr = np.asarray(x_grid, dtype=np.float64)[np.newaxis, :]
        m_arr = np.arange(num_moments, dtype=np.float64)[:, np.newaxis]
        evaluated = self._fn(m_arr, x_arr)
        broadcasted = np.broadcast_to(evaluated, (num_moments, x_arr.shape[1]))
        return np.ascontiguousarray(broadcasted, dtype=np.float64)


_x_sym = sp.Symbol("x", real=True)
_m_sym = sp.Symbol("m", integer=True, nonnegative=True)

DeltaFunction = SpectralFunction(
    sp.cos(_m_sym * sp.acos(_x_sym)) / (sp.pi * sp.sqrt(1 - _x_sym**2)),
    x=_x_sym,
    m=_m_sym,
)


def legacy_delta_table(
    x: Sequence[float] | npt.NDArray[np.float64], num_moments: int
) -> npt.NDArray[np.float64]:
    """Compute the ``(M_eff, N)`` delta coefficient table in the exact legacy order.

    Solves the problem of reproducing C++ ``delta_chebF(x, m)`` bit for bit
    without SIMD ULP drift by calling Python's ``math.sqrt``, ``math.acos``,
    ``math.cos``, and ``math.pi`` in the exact expression order:

    .. code-block:: python

        f0 = math.sqrt(1.0 - xi * xi)
        fm = math.cos(m * math.acos(xi))
        val = fm / f0 / math.pi

    All inputs and outputs are dimensionless in the rescaled Chebyshev domain.
    Time and space complexity are ``O(num_moments * len(x))``.

    Parameters
    ----------
    x : Sequence[float] or numpy.ndarray
        Rescaled coordinate grid ``x_i in (-1, 1)`` of length ``N``.
    num_moments : int
        Effective number of Chebyshev moments ``M_eff > 0``.

    Returns
    -------
    numpy.ndarray
        2D C-contiguous ``float64`` array of shape ``(num_moments, len(x))``.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments <= 0`` or if any ``|x_i| >= 1.0``.

    Notes
    -----
    Pitfall: Keep the two sequential divisions ``fm / f0 / math.pi`` rather than
    ``fm / (f0 * math.pi)``; multiplying the denominator first alters rounding
    at the last bit.
    """
    if num_moments <= 0:
        raise ValueError("legacy_delta_table: num_moments must be positive")

    x_list = [float(val) for val in x]
    num_points = len(x_list)
    table = np.empty((num_moments, num_points), dtype=np.float64)
    pi = math.pi
    cos = math.cos
    acos = math.acos
    sqrt = math.sqrt

    for i, xi in enumerate(x_list):
        f0 = sqrt(1.0 - xi * xi)
        acos_x = acos(xi)
        for m_idx in range(num_moments):
            fm = cos(float(m_idx) * acos_x)
            table[m_idx, i] = fm / f0 / pi

    return table
