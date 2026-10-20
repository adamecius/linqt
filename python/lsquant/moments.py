"""Bytecode plan constructors, spectral window helpers, and execution routines.

Defines the rank-1 density-of-states plan (:func:`dos_plan`) and generalised
rank-2 Kubo correlation plan (:func:`kubo_plan`) as data over the
five-instruction C++ engine grammar, together with
:func:`window_from_gershgorin` and :func:`run_plan`.
"""

from __future__ import annotations

from typing import Sequence
from . import _core


def dos_plan(num_moments: int, num_states: int) -> _core.Plan:
    """Build the rank-1 density-of-states bytecode :class:`~lsquant._core.Plan`.

    Solves the problem of expressing ``Tr[1 * delta(E - H)]`` as a ket-only
    Chebyshev expansion ``Expand(M)`` followed by ``Trace(w)`` with the legacy
    per-state KPM weights ``w[0] = 1 / num_states`` and
    ``w[m] = 2 / num_states`` for ``m > 0``. Complexity is ``O(num_moments)``.

    Parameters
    ----------
    num_moments : int
        Number of Chebyshev moments ``M > 0``.
    num_states : int
        Number of trace vectors ``N_states > 0`` used for weight normalisation.

    Returns
    -------
    lsquant._core.Plan
        Compiled bytecode plan with ``word = "1 delta"``.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments <= 0`` or ``num_states <= 0``.

    Notes
    -----
    Pitfall: Because the ``(2 - delta_{m,0}) / num_states`` weights are folded
    into :class:`~lsquant._core.Trace` per state, do not apply KPM weights a
    second time after running this plan.
    """
    if num_moments <= 0:
        raise ValueError("dos_plan: num_moments must be positive")
    if num_states <= 0:
        raise ValueError("dos_plan: num_states must be positive")

    weights: list[complex] = []
    base = 2.0 / float(num_states)
    for m in range(num_moments):
        scal = base * 0.5 if m == 0 else base
        weights.append(complex(scal, 0.0))

    code = [
        _core.Expand(num_moments),
        _core.Trace(weights),
    ]
    return _core.Plan(code, "1 delta")


def kubo_plan(
    num_moments: int,
    left_index: int = 0,
    right_index: int = 0,
) -> _core.Plan:
    """Build the rank-2 Kubo correlation bytecode :class:`~lsquant._core.Plan`.

    Solves the problem of expressing the two-operator correlation trace
    ``mu_{m,n} = Tr[O_left T_m(H~) O_right T_n(H~)]`` in the exact ket/bra
    evaluation order of legacy ``CorrelationExpansionMoments(OPL, OPR)``:

    - ket: ``Expand(M)``, ``ApplyOperator(right_index)``
      (computes ``O_right T_n |r>``);
    - ``MirrorToBra()``;
    - bra: ``ApplyOperator(left_index)``, ``Expand(M)``
      (computes ``T_m O_left |r>``);
    - ``Trace([1])`` (accumulates raw inner products across states).

    Complexity is ``O(1)``.

    Parameters
    ----------
    num_moments : int
        Number of Chebyshev moments ``M > 0`` along each of the two axes.
        (If called as ``kubo_plan(left_index, right_index, num_moments)`` with
        small operator indices ``0, 1`` first and ``M >= 2`` third, the order
        is normalised automatically.)
    left_index : int, optional
        Index of the bra-side operator ``O_left`` in the engine's operator list
        (default ``0``).
    right_index : int, optional
        Index of the ket-side operator ``O_right`` in the engine's operator list
        (default ``0``).

    Returns
    -------
    lsquant._core.Plan
        Compiled bytecode plan.

    Warnings
    --------
    Raises ``ValueError`` if ``num_moments <= 0`` or operator indices are
    negative.

    Notes
    -----
    Pitfall: In legacy ``inline_compute-kpm-nonEqOp <label> <OP_R> <OP_L> <M>``,
     the first CLI operator seeds the bra recursion (``left_index``) and the
    second is applied to the ket batch (``right_index``).
    """
    if num_moments <= 1 and right_index >= 2:
        num_moments, left_index, right_index = right_index, num_moments, left_index

    if num_moments <= 0:
        raise ValueError("kubo_plan: num_moments must be positive")
    if left_index < 0 or right_index < 0:
        raise ValueError("kubo_plan: operator indices must be non-negative")

    code = [
        _core.Expand(num_moments),
        _core.ApplyOperator(right_index),
        _core.MirrorToBra(),
        _core.ApplyOperator(left_index),
        _core.Expand(num_moments),
        _core.Trace([complex(1.0, 0.0)]),
    ]
    word = "VX delta VX delta" if left_index == right_index else "VX delta VY delta"
    return _core.Plan(code, word)


def window_from_gershgorin(
    hamiltonian: _core.SparseMatrix, pad: float = 0.10
) -> _core.SpectralWindow:
    """Compute the symmetric Gershgorin :class:`~lsquant._core.SpectralWindow` for ``hamiltonian``.

    Solves the problem of enclosing the spectrum of ``hamiltonian`` when no
    explicit ``BOUNDS`` or ``.desc`` sidecar is supplied, reproducing the legacy
    ``utility::SpectralBounds`` fallback ``[-rho * (1 + pad), +rho * (1 + pad)]``
    where ``rho = max_i sum_j |H_{ij}|`` is the infinity norm of ``hamiltonian``.
    Complexity is ``O(dim + nnz)``.

    Parameters
    ----------
    hamiltonian : lsquant._core.SparseMatrix
        Compressed sparse Hamiltonian matrix.
    pad : float, optional
        Relative safety padding fraction (default ``0.10``).

    Returns
    -------
    lsquant._core.SpectralWindow
        Symmetric spectral window ``[-rho * (1 + pad), +rho * (1 + pad)]``.

    Warnings
    --------
    Raises ``ValueError`` if ``pad < 0``.

    Notes
    -----
    Pitfall: Ensure ``hamiltonian`` has undergone :meth:`hermitian_fold` (as
    :func:`~lsquant._core.ReadCsrText` does) before computing Gershgorin row
    sums so that the resulting ``BandWidth`` matches legacy files to the last
    bit.
    """
    bounds = _core.GershgorinBounds(hamiltonian, float(pad))
    return _core.SpectralWindow(bounds.lower, bounds.upper)


def run_plan(
    plan: _core.Plan,
    hamiltonian: _core.SparseMatrix,
    window: _core.SpectralWindow,
    states: _core.StateGenerator,
    operators: Sequence[_core.LinearOperator] = (),
    batch_size: int | None = None,
    symmetrise: bool = True,
) -> _core.ChebyshevMoments:
    """Execute a bytecode :class:`~lsquant._core.Plan` and populate tensor metadata.

    Rescales ``hamiltonian`` to ``window``, runs ``plan`` through
    :class:`~lsquant._core.MomentEngine`, and records ``window``,
    ``kpm_weighted``, and ``symmetrised`` on the returned
    :class:`~lsquant._core.ChebyshevMoments`. For rank-2 plans with unit trace
    weight, applies :func:`lsquant._core.ApplyKpmWeights`.

    Parameters
    ----------
    plan : lsquant._core.Plan
        Compiled bytecode instruction sequence.
    hamiltonian : lsquant._core.SparseMatrix
        Physical tight-binding Hamiltonian in CSR format.
    window : lsquant._core.SpectralWindow
        Spectral window ``[lower, upper]`` enclosing the spectrum of ``hamiltonian``.
    states : lsquant._core.StateGenerator
        Source-state generator for trace evaluation.
    operators : Sequence[lsquant._core.LinearOperator], optional
        Positional list of static operators referenced by ``ApplyOperator``.
    batch_size : int or None, optional
        Chebyshev batch size (defaults to the first ``Expand`` moment count).
    symmetrise : bool, optional
        Whether to Hermitian-symmetrise rank-2 moment tensors (default ``True``).

    Returns
    -------
    lsquant._core.ChebyshevMoments
        Accumulated and weighted Chebyshev moment tensor.

    Warnings
    --------
    Raises ``ValueError`` or ``RuntimeError`` if dimensions mismatch or ``plan``
    is invalid.

    Notes
    -----
    Pitfall: ``states`` is automatically reset by ``MomentEngine.run`` before
    traversing the state sequence.
    """
    rescaled_h = _core.RescaledOperator(hamiltonian, window)
    if batch_size is None:
        batch_size = 1
        for instr in plan.code:
            if isinstance(instr, _core.Expand):
                batch_size = instr.num_moments
                break

    engine = _core.MomentEngine(rescaled_h, list(operators), states, batch_size)
    mu = engine.run(plan)
    mu.window = window
    mu.system_size = hamiltonian.dimension
    mu.number_of_states = int(states.number_of_states)
    mu.state_label = states.label
    mu.word = plan.word

    if mu.rank == 1:
        mu.kpm_weighted = True
        mu.symmetrised = False
    elif mu.rank == 2 and not mu.kpm_weighted:
        _core.ApplyKpmWeights(mu, int(states.number_of_states), symmetrise)

    return mu
