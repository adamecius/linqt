"""Basis adapters (RealBasis, BlochBasis, EigenBasis) for the Chebyshev engine.

All three bases expose the same interface over a :class:`~lsquant.models.wannier.WannierModel`
and supercell ``(N_1, N_2, N_3)``:

- ``.hamiltonian() -> LinearOperator``
- ``.operator(name) -> LinearOperator``
- ``.states(trace_mode, num_states=None, seed=12345) -> StateGenerator``
- ``.window(pad=0.10) -> SpectralWindow``

The C++ Chebyshev engine (:class:`~lsquant._core.MomentEngine`) receives only
:class:`~lsquant._core.LinearOperator` and :class:`~lsquant._core.StateGenerator`
objects and is completely agnostic of which basis adapter produced them.
"""

from __future__ import annotations

import json
import math
import os
from typing import Sequence
import numpy as np
import numpy.typing as npt

from . import _core
from .models.wannier import WannierModel
from .moments import window_from_gershgorin
from .spectral_functions import GaussianDelta


def make_anderson_disorder(
    dimension: int, disorder_w: float, seed: int = 12345
) -> npt.NDArray[np.complex128]:
    """Generate deterministic real-space box disorder in ``[-W/2, +W/2]``.

    Parameters
    ----------
    dimension : int
        Total Hilbert-space dimension ``N_k * W``.
    disorder_w : float
        Box disorder width ``W >= 0`` in eV.
    seed : int, optional
        Deterministic PRNG seed (default ``12345``).

    Returns
    -------
    numpy.ndarray
        1D ``complex128`` array of length ``dimension`` with real entries drawn
        uniformly from ``[-W/2, +W/2]``.
    """
    if disorder_w < 0.0:
        raise ValueError("make_anderson_disorder: disorder_w must be non-negative")
    if disorder_w == 0.0:
        return np.zeros(dimension, dtype=np.complex128)
    rng = np.random.default_rng(int(seed))
    vals = rng.uniform(-0.5 * float(disorder_w), 0.5 * float(disorder_w), size=dimension)
    return np.ascontiguousarray(vals, dtype=np.complex128)


class RealBasis:
    """Real-space supercell basis adapter.

    Reads ``<seed>.HAM.CSR`` and ``<seed>.VX.CSR`` (and ``VY``, ``VZ``) written
    by ``wannier2sparse`` when they exist in ``csr_dir`` or ``model.model_dir``
    for the requested supercell ``N_1 * N_2 * N_3 * W``, and builds them
    identically in the ``wannier2sparse`` cell numbering
    ``iR = n_1 + N_1 * n_2 + N_1 * N_2 * n_3`` otherwise. When ``disorder_w > 0``
    (or an explicit ``disorder`` vector is given), ``H`` is the matrix-free
    composition ``OperatorSum(SparseMatrixOperator(H_csr), DiagonalOperator(v_real))``.
    """

    def __init__(
        self,
        model: WannierModel,
        supercell: Sequence[int],
        *,
        disorder_w: float = 0.0,
        seed: int = 12345,
        disorder_seed: int | None = None,
        disorder: npt.NDArray[np.complex128] | npt.NDArray[np.float64] | None = None,
        csr_dir: str | None = None,
        use_w2s_csr: bool = True,
    ) -> None:
        if len(supercell) != 3:
            raise ValueError("RealBasis: supercell must have 3 dimensions (N1, N2, N3)")
        self.model = model
        self.supercell = (int(supercell[0]), int(supercell[1]), int(supercell[2]))
        self.num_cells = self.supercell[0] * self.supercell[1] * self.supercell[2]
        self.num_wann = model.num_wann
        self.dimension = self.num_cells * self.num_wann
        self.disorder_w = float(disorder_w)
        self.seed = int(seed)
        self.disorder_seed = self.seed if disorder_seed is None else int(disorder_seed)
        self._csr_dir = csr_dir
        self._use_w2s_csr = bool(use_w2s_csr)

        loaded_ham = self._try_load_w2s_csr("HAM")
        self._h_csr: _core.SparseMatrix = (
            loaded_ham
            if loaded_ham is not None
            else model.hamiltonian_csr(self.supercell)
        )
        self._h_op_clean = _core.SparseMatrixOperator(self._h_csr)

        if disorder is not None:
            v_arr = np.ascontiguousarray(disorder, dtype=np.complex128).ravel()
            if v_arr.size != self.dimension:
                raise ValueError("RealBasis: disorder array length must equal dimension")
            self._disorder_vec: npt.NDArray[np.complex128] | None = v_arr
            self._v_op: _core.DiagonalOperator | None = _core.DiagonalOperator(v_arr)
            self._h_op: _core.LinearOperator = _core.OperatorSum(
                self._h_op_clean, self._v_op
            )
        elif self.disorder_w > 0.0:
            v_arr = make_anderson_disorder(
                self.dimension, self.disorder_w, self.disorder_seed
            )
            self._disorder_vec = v_arr
            self._v_op = _core.DiagonalOperator(v_arr)
            self._h_op = _core.OperatorSum(self._h_op_clean, self._v_op)
        else:
            self._disorder_vec = None
            self._v_op = None
            self._h_op = self._h_op_clean

        self._vel_csr_cache: dict[str, _core.SparseMatrix] = {}
        self._op_cache: dict[str, _core.LinearOperator] = {}

    def _try_load_w2s_csr(self, op_name: str) -> _core.SparseMatrix | None:
        """Load ``<seed>.<op_name>.CSR`` written by ``wannier2sparse`` if present and matching ``self.supercell``."""
        if not self._use_w2s_csr:
            return None
        seed = self.model.seed_name
        search_dirs: list[str] = []
        if self._csr_dir is not None:
            search_dirs.append(self._csr_dir)
        if self.model.model_dir is not None and self.model.model_dir not in search_dirs:
            search_dirs.append(self.model.model_dir)
            ops_sub = os.path.join(self.model.model_dir, "operators")
            if os.path.isdir(ops_sub) and ops_sub not in search_dirs:
                search_dirs.append(ops_sub)

        for d in search_dirs:
            csr_path = os.path.join(d, f"{seed}.{op_name}.CSR")
            if not os.path.isfile(csr_path):
                continue
            out_path = os.path.join(d, f"{seed}.out")
            if os.path.isfile(out_path):
                try:
                    with open(out_path, encoding="utf-8") as fh:
                        meta = json.load(fh)
                    sc = meta.get("supercell")
                    if isinstance(sc, list) and len(sc) == 3:
                        if (int(sc[0]), int(sc[1]), int(sc[2])) != self.supercell:
                            continue
                except Exception:  # noqa: BLE001
                    pass
            try:
                mat = _core.ReadCsrText(csr_path)
            except Exception:  # noqa: BLE001
                continue
            if int(mat.dimension) == self.dimension:
                return mat
        return None

    @property
    def hamiltonian_csr(self) -> _core.SparseMatrix:
        """Return the clean real-space supercell Hamiltonian :class:`~lsquant._core.SparseMatrix`."""
        return self._h_csr

    @property
    def disorder_values(self) -> npt.NDArray[np.complex128] | None:
        """Return the real-space diagonal disorder array if present."""
        return self._disorder_vec

    def hamiltonian(self) -> _core.LinearOperator:
        """Return the Hamiltonian :class:`~lsquant._core.LinearOperator` in real space."""
        return self._h_op

    def velocity_csr(self, name: str) -> _core.SparseMatrix:
        """Return the real-space supercell velocity :class:`~lsquant._core.SparseMatrix`."""
        key = name.strip().upper()
        if key not in ("VX", "VY", "VZ"):
            raise ValueError(f"RealBasis.velocity_csr: unsupported velocity '{name}'")
        if key not in self._vel_csr_cache:
            loaded_vel = self._try_load_w2s_csr(key)
            self._vel_csr_cache[key] = (
                loaded_vel
                if loaded_vel is not None
                else self.model.velocity_csr(self.supercell, key)
            )
        return self._vel_csr_cache[key]

    def operator(self, name: str) -> _core.LinearOperator:
        """Return the named static operator as a :class:`~lsquant._core.LinearOperator`."""
        key = name.strip().upper()
        if key in ("1", "I", "ID", "IDENTITY"):
            if "1" not in self._op_cache:
                self._op_cache["1"] = _core.IdentityOperator(self.dimension)
            return self._op_cache["1"]
        if key in ("H", "HAM"):
            return self._h_op
        if key in ("VX", "VY", "VZ"):
            if key not in self._op_cache:
                mat = self.velocity_csr(key)
                self._op_cache[key] = _core.SparseMatrixOperator(mat)
            return self._op_cache[key]
        raise ValueError(f"RealBasis.operator: unknown operator '{name}'")

    def states(
        self,
        trace_mode: str = "exact",
        num_states: int | None = None,
        seed: int | None = None,
    ) -> _core.StateGenerator:
        """Return a :class:`~lsquant._core.StateGenerator` for ``trace_mode``."""
        mode = trace_mode.strip().lower()
        rng_seed = self.seed if seed is None else int(seed)
        if mode == "exact":
            n_st = self.dimension if num_states is None else int(num_states)
            if n_st <= 0 or n_st > self.dimension:
                raise ValueError(
                    "RealBasis.states: num_states must lie in 1..dimension for exact trace"
                )
            return _core.LocalStates(self.dimension, list(range(n_st)))
        if mode == "random":
            n_st = 1 if num_states is None else int(num_states)
            if n_st <= 0:
                raise ValueError("RealBasis.states: num_states must be positive")
            return _core.RandomPhaseStates(self.dimension, n_st, rng_seed)
        raise ValueError(
            f"RealBasis.states: unsupported trace_mode '{trace_mode}' (expected 'exact' or 'random')"
        )

    def window(self, pad: float = 0.10) -> _core.SpectralWindow:
        """Return the Gershgorin :class:`~lsquant._core.SpectralWindow` for this Hamiltonian.

        Kept on Gershgorin only as the stage-1-compatible real-space default;
        on Wannier models with long-range hoppings, Gershgorin bounds are
        several times wider than the true spectrum and waste Chebyshev moments,
        so :class:`BlochBasis` and :class:`EigenBasis` use the block eigenvalue
        range instead.
        """
        win = window_from_gershgorin(self._h_csr, pad=float(pad))
        if self.disorder_w > 0.0:
            extra = 0.5 * self.disorder_w * (1.0 + float(pad))
            return _core.SpectralWindow(win.lower - extra, win.upper + extra)
        return win


class BlochBasis:
    """Bloch k-space basis adapter with optional FFT-sandwiched real-space disorder.

    Builds ``H_k = BlockDiagonalOperator(H(k))`` and ``v_a(k) = BlockDiagonalOperator(dH_a(k))``
    over the ``N_k = N_1 * N_2 * N_3`` supercell k-points, and owns the batched
    FFTW :class:`~lsquant._core.BlochTransform` ``B`` (`real -> Bloch`). When
    ``disorder_w > 0`` (or ``disorder`` is supplied), ``H`` is the composite
    ``OperatorSum(H_k, SandwichOperator(B, DiagonalOperator(v_real)))``.
    """

    def __init__(
        self,
        model: WannierModel,
        supercell: Sequence[int],
        *,
        disorder_w: float = 0.0,
        seed: int = 12345,
        disorder_seed: int | None = None,
        disorder: npt.NDArray[np.complex128] | npt.NDArray[np.float64] | None = None,
    ) -> None:
        if len(supercell) != 3:
            raise ValueError("BlochBasis: supercell must have 3 dimensions (N1, N2, N3)")
        self.model = model
        self.supercell = (int(supercell[0]), int(supercell[1]), int(supercell[2]))
        self.num_cells = self.supercell[0] * self.supercell[1] * self.supercell[2]
        self.num_wann = model.num_wann
        self.dimension = self.num_cells * self.num_wann
        self.disorder_w = float(disorder_w)
        self.seed = int(seed)
        self.disorder_seed = self.seed if disorder_seed is None else int(disorder_seed)

        self.k_points = model.k_mesh(self.supercell)
        self.phases = model.bloch_phases(self.k_points)
        self.transform = _core.BlochTransform(
            self.supercell[0],
            self.supercell[1],
            self.supercell[2],
            self.num_wann,
            self.phases.ravel(),
            _core.BlochDirection.RealToBloch,
        )

        self.hk_array = model.Hk(self.k_points)
        self._hk_op = _core.BlockDiagonalOperator(self.hk_array)

        if disorder is not None:
            v_arr = np.ascontiguousarray(disorder, dtype=np.complex128).ravel()
            if v_arr.size != self.dimension:
                raise ValueError("BlochBasis: disorder array length must equal dimension")
            self._disorder_vec: npt.NDArray[np.complex128] | None = v_arr
            self._v_op: _core.DiagonalOperator | None = _core.DiagonalOperator(v_arr)
            self._sandwich_v: _core.SandwichOperator | None = _core.SandwichOperator(
                self.transform, self._v_op
            )
            self._h_op: _core.LinearOperator = _core.OperatorSum(
                self._hk_op, self._sandwich_v
            )
        elif self.disorder_w > 0.0:
            v_arr = make_anderson_disorder(
                self.dimension, self.disorder_w, self.disorder_seed
            )
            self._disorder_vec = v_arr
            self._v_op = _core.DiagonalOperator(v_arr)
            self._sandwich_v = _core.SandwichOperator(self.transform, self._v_op)
            self._h_op = _core.OperatorSum(self._hk_op, self._sandwich_v)
        else:
            self._disorder_vec = None
            self._v_op = None
            self._sandwich_v = None
            self._h_op = self._hk_op

        self._op_cache: dict[str, _core.LinearOperator] = {}
        self._eig_bounds: tuple[float, float] | None = None

    @property
    def disorder_values(self) -> npt.NDArray[np.complex128] | None:
        """Return the real-space diagonal disorder array if present."""
        return self._disorder_vec

    def hamiltonian(self) -> _core.LinearOperator:
        """Return the Hamiltonian :class:`~lsquant._core.LinearOperator` in Bloch space."""
        return self._h_op

    def operator(self, name: str) -> _core.LinearOperator:
        """Return the named static operator in Bloch space."""
        key = name.strip().upper()
        if key in ("1", "I", "ID", "IDENTITY"):
            if "1" not in self._op_cache:
                self._op_cache["1"] = _core.IdentityOperator(self.dimension)
            return self._op_cache["1"]
        if key in ("H", "HAM"):
            return self._h_op
        if key in ("VX", "VY", "VZ"):
            if key not in self._op_cache:
                vel_blocks = self.model.dHk(self.k_points, key)
                self._op_cache[key] = _core.BlockDiagonalOperator(vel_blocks)
            return self._op_cache[key]
        raise ValueError(f"BlochBasis.operator: unknown operator '{name}'")

    def states(
        self,
        trace_mode: str = "exact",
        num_states: int | None = None,
        seed: int | None = None,
    ) -> _core.StateGenerator:
        """Return a :class:`~lsquant._core.StateGenerator` for ``trace_mode``.

        For ``trace_mode == "exact"`` at ``disorder_w == 0``, returns
        :class:`~lsquant._core.BlockLocalStates` ``(N_k, W)``, which evaluates
        the exact full-mesh trace of any block-diagonal word in ``W`` states.
        For ``trace_mode == "random"``, applies ``B`` to
        :class:`~lsquant._core.RandomPhaseStates` via
        :class:`~lsquant._core.TransformedStates` so that real-space and
        Bloch-space stochastic traces use the exact same physical state vectors.
        """
        mode = trace_mode.strip().lower()
        rng_seed = self.seed if seed is None else int(seed)
        if mode == "exact":
            if self._disorder_vec is None and (
                num_states is None or int(num_states) == self.num_wann
            ):
                return _core.BlockLocalStates(self.num_cells, self.num_wann)
            n_st = self.dimension if num_states is None else int(num_states)
            if n_st <= 0 or n_st > self.dimension:
                raise ValueError(
                    "BlochBasis.states: num_states must lie in 1..dimension for full exact trace"
                )
            return _core.LocalStates(self.dimension, list(range(n_st)))
        if mode == "random":
            n_st = 1 if num_states is None else int(num_states)
            if n_st <= 0:
                raise ValueError("BlochBasis.states: num_states must be positive")
            inner = _core.RandomPhaseStates(self.dimension, n_st, rng_seed)
            return _core.TransformedStates(self.transform, inner)
        raise ValueError(
            f"BlochBasis.states: unsupported trace_mode '{trace_mode}' (expected 'exact' or 'random')"
        )

    def window(self, pad: float = 0.10) -> _core.SpectralWindow:
        """Return the :class:`~lsquant._core.SpectralWindow` enclosing the block eigenvalues.

        Diagonalises the ``N_k`` blocks ``H(k)`` once and caches ``(emin, emax)``,
        returning ``[emin - pad * width, emax + pad * width]`` with
        ``width = emax - emin`` identical to :meth:`EigenBasis.window`.
        On Wannier models with long-range hoppings, Gershgorin bounds are
        several times wider than the true spectrum and waste Chebyshev moments,
        so ``BlochBasis`` uses the block eigenvalues rather than Gershgorin.
        """
        if pad < 0.0:
            raise ValueError("BlochBasis.window: pad must be non-negative")
        if self._eig_bounds is None:
            eigvals, _ = np.linalg.eigh(self.hk_array)
            self._eig_bounds = (float(np.min(eigvals)), float(np.max(eigvals)))
        emin, emax = self._eig_bounds
        if self.disorder_w > 0.0:
            emin -= 0.5 * self.disorder_w
            emax += 0.5 * self.disorder_w
        width = emax - emin
        if width <= 0.0:
            emin -= 1.0
            emax += 1.0
            width = 2.0
        margin = float(pad) * width
        return _core.SpectralWindow(emin - margin, emax + margin)


class EigenBasis:
    """Per-k exact eigenbasis adapter (clean systems, ``W = 0``).

    Diagonalises each ``W x W`` block ``H(k) = U(k) diag(eps(k)) U(k)^+`` via
    :func:`numpy.linalg.eigh`, builds ``H = DiagonalOperator(eps)`` and rotated
    operators ``O_eig(k) = U(k)^+ O(k) U(k)`` as
    :class:`~lsquant._core.BlockDiagonalOperator`, and uses
    :class:`~lsquant._core.BlockLocalStates` ``(N_k, W)`` so the standard
    :class:`~lsquant._core.MomentEngine` evaluates exact full-mesh moments in
    ``W`` state evaluations.
    """

    def __init__(
        self,
        model: WannierModel,
        supercell: Sequence[int],
        *,
        disorder_w: float = 0.0,
        seed: int = 12345,
    ) -> None:
        if float(disorder_w) != 0.0:
            raise ValueError(
                "EigenBasis: non-zero disorder_w is not supported in the per-k eigenbasis"
            )
        if len(supercell) != 3:
            raise ValueError("EigenBasis: supercell must have 3 dimensions (N1, N2, N3)")
        self.model = model
        self.supercell = (int(supercell[0]), int(supercell[1]), int(supercell[2]))
        self.num_cells = self.supercell[0] * self.supercell[1] * self.supercell[2]
        self.num_wann = model.num_wann
        self.dimension = self.num_cells * self.num_wann
        self.disorder_w = 0.0
        self.seed = int(seed)

        self.k_points = model.k_mesh(self.supercell)
        self.hk_array = model.Hk(self.k_points)

        # Batched eigh over all N_k blocks of shape (N_k, W, W)
        eigvals, eigvecs = np.linalg.eigh(self.hk_array)
        self.eigenvalues: npt.NDArray[np.float64] = np.ascontiguousarray(
            eigvals, dtype=np.float64
        )
        self.eigenvectors: npt.NDArray[np.complex128] = np.ascontiguousarray(
            eigvecs, dtype=np.complex128
        )
        self.eigenvectors_dagger: npt.NDArray[np.complex128] = np.ascontiguousarray(
            np.conj(np.swapaxes(self.eigenvectors, -1, -2)), dtype=np.complex128
        )

        self._h_op = _core.DiagonalOperator(
            np.ascontiguousarray(self.eigenvalues.ravel(), dtype=np.complex128)
        )
        self.u_op = _core.BlockDiagonalOperator(self.eigenvectors)
        self.u_dagger_op = _core.BlockDiagonalOperator(self.eigenvectors_dagger)
        self._op_cache: dict[str, _core.LinearOperator] = {}

    def hamiltonian(self) -> _core.LinearOperator:
        """Return ``DiagonalOperator(eps)``."""
        return self._h_op

    def rotate_blocks(
        self, bloch_blocks: npt.NDArray[np.complex128]
    ) -> npt.NDArray[np.complex128]:
        """Rotate ``(N_k, W, W)`` Bloch blocks into the eigenbasis: ``U(k)^+ O(k) U(k)``."""
        rotated = self.eigenvectors_dagger @ bloch_blocks @ self.eigenvectors
        return np.ascontiguousarray(rotated, dtype=np.complex128)

    def operator(self, name: str) -> _core.LinearOperator:
        """Return the named static operator rotated into the eigenbasis."""
        key = name.strip().upper()
        if key in ("1", "I", "ID", "IDENTITY"):
            if "1" not in self._op_cache:
                self._op_cache["1"] = _core.IdentityOperator(self.dimension)
            return self._op_cache["1"]
        if key in ("H", "HAM"):
            return self._h_op
        if key in ("VX", "VY", "VZ"):
            if key not in self._op_cache:
                vel_bloch = self.model.dHk(self.k_points, key)
                vel_eig = self.rotate_blocks(vel_bloch)
                self._op_cache[key] = _core.BlockDiagonalOperator(vel_eig)
            return self._op_cache[key]
        raise ValueError(f"EigenBasis.operator: unknown operator '{name}'")

    def states(
        self,
        trace_mode: str = "exact",
        num_states: int | None = None,
        seed: int | None = None,
    ) -> _core.StateGenerator:
        """Return a :class:`~lsquant._core.StateGenerator` for ``trace_mode``.

        Defaults to :class:`~lsquant._core.BlockLocalStates` ``(N_k, W)``, which
        evaluates the exact eigenbasis trace in ``W`` states.
        """
        mode = trace_mode.strip().lower()
        rng_seed = self.seed if seed is None else int(seed)
        if mode == "exact":
            if num_states is None or int(num_states) == self.num_wann:
                return _core.BlockLocalStates(self.num_cells, self.num_wann)
            n_st = int(num_states)
            if n_st <= 0 or n_st > self.dimension:
                raise ValueError(
                    "EigenBasis.states: num_states must lie in 1..dimension"
                )
            return _core.LocalStates(self.dimension, list(range(n_st)))
        if mode == "random":
            n_st = 1 if num_states is None else int(num_states)
            if n_st <= 0:
                raise ValueError("EigenBasis.states: num_states must be positive")
            phases = self.model.bloch_phases(self.k_points)
            b_op = _core.BlochTransform(
                self.supercell[0],
                self.supercell[1],
                self.supercell[2],
                self.num_wann,
                phases.ravel(),
                _core.BlochDirection.RealToBloch,
            )
            u_dag_b = _core.OperatorProduct(self.u_dagger_op, b_op)
            inner = _core.RandomPhaseStates(self.dimension, n_st, rng_seed)
            return _core.TransformedStates(u_dag_b, inner)
        raise ValueError(
            f"EigenBasis.states: unsupported trace_mode '{trace_mode}' (expected 'exact' or 'random')"
        )

    def window(self, pad: float = 0.10) -> _core.SpectralWindow:
        """Return the :class:`~lsquant._core.SpectralWindow` enclosing the block eigenvalues.

        Computes ``[emin - pad * width, emax + pad * width]`` from the minimum and
        maximum eigenvalues ``(emin, emax)`` of all ``N_k`` blocks ``H(k)`` with
        ``width = emax - emin``. On Wannier models with long-range hoppings,
        Gershgorin bounds are several times wider than the true spectrum and
        waste Chebyshev moments, so ``EigenBasis`` and :class:`BlochBasis`
        take their window from the block eigenvalues rather than from Gershgorin.
        """
        if pad < 0.0:
            raise ValueError("EigenBasis.window: pad must be non-negative")
        emin = float(np.min(self.eigenvalues))
        emax = float(np.max(self.eigenvalues))
        width = emax - emin
        if width <= 0.0:
            emin -= 1.0
            emax += 1.0
            width = 2.0
        margin = float(pad) * width
        return _core.SpectralWindow(emin - margin, emax + margin)

    def exact_dos(
        self,
        energies: Sequence[float] | npt.NDArray[np.float64],
        eta: float | Sequence[float] | npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate the exact per-cell Gaussian density of states on the block eigenvalues.

        Computes ``(1 / N_k) sum_{k, n} g(E - eps_n(k); eta)`` where the
        Gaussian ``g(x; eta) = exp(-x^2 / (2 * eta^2)) / (sqrt(2 * pi) * eta)``
        is evaluated via the stage-2 SymPy
        :data:`~lsquant.spectral_functions.GaussianDelta`
        :class:`~lsquant.spectral_functions.SpectralFunction`.

        Parameters
        ----------
        energies : Sequence[float] or numpy.ndarray
            1D energy grid in eV.
        eta : float or numpy.ndarray
            Gaussian broadening width in eV (``eta > 0``).

        Returns
        -------
        numpy.ndarray
            1D ``float64`` array of length ``len(energies)`` giving the density
            of states per unit cell in ``1 / eV``.
        """
        e_arr = np.asarray(energies, dtype=np.float64).ravel()
        eta_arr = np.asarray(eta, dtype=np.float64)
        if np.any(eta_arr <= 0.0):
            raise ValueError("EigenBasis.exact_dos: eta must be positive")
        if eta_arr.ndim == 1:
            eta_arr = eta_arr[:, None, None]
        diff = e_arr[:, None, None] - self.eigenvalues[None, :, :]  # (N_E, N_k, W)
        g_vals = np.asarray(GaussianDelta.evaluate(diff, eta_arr), dtype=np.float64)
        dos_per_cell = np.sum(g_vals, axis=(1, 2)) / float(self.num_cells)
        return np.ascontiguousarray(dos_per_cell, dtype=np.float64)

    def exact_sigma_xx(
        self,
        energies: Sequence[float] | npt.NDArray[np.float64],
        eta: float | Sequence[float] | npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate the exact Gaussian Kubo-Greenwood conductivity ``sigma_xx(E)``.

        Computes ``pi * sum_{k, n, m} |v_{x, nm}(k)|^2 g(E - eps_n(k); eta) g(E - eps_m(k); eta)``
        with the same normalisation as stage-3 :func:`~lsquant.reconstruct.kubo_greenwood`
        so that the two curves are directly comparable (multiply by
        ``2 * pi / (A_cell * N_k)`` to convert to ``hr_exactdiag.py``'s
        ``e^2 / h`` per-area convention).

        Parameters
        ----------
        energies : Sequence[float] or numpy.ndarray
            1D energy grid in eV.
        eta : float or numpy.ndarray
            Gaussian broadening width in eV (``eta > 0``).

        Returns
        -------
        numpy.ndarray
            1D ``float64`` array of length ``len(energies)``.
        """
        e_arr = np.asarray(energies, dtype=np.float64).ravel()
        eta_arr = np.asarray(eta, dtype=np.float64)
        if np.any(eta_arr <= 0.0):
            raise ValueError("EigenBasis.exact_sigma_xx: eta must be positive")
        if eta_arr.ndim == 1:
            eta_arr = eta_arr[:, None, None]
        vx_eig = self.rotate_blocks(self.model.dHk(self.k_points, "VX"))
        vx2 = np.abs(vx_eig) ** 2  # (N_k, W, W)
        diff = e_arr[:, None, None] - self.eigenvalues[None, :, :]  # (N_E, N_k, W)
        g_vals = np.asarray(GaussianDelta.evaluate(diff, eta_arr), dtype=np.float64)
        sigma_ext = math.pi * np.einsum(
            "ekn,knm,ekm->e", g_vals, vx2, g_vals, optimize=True
        )
        return np.ascontiguousarray(sigma_ext, dtype=np.float64)


def make_basis(
    basis_name: str,
    model: WannierModel,
    supercell: Sequence[int],
    *,
    disorder_w: float = 0.0,
    seed: int = 12345,
    disorder_seed: int | None = None,
    disorder: npt.NDArray[np.complex128] | npt.NDArray[np.float64] | None = None,
    csr_dir: str | None = None,
) -> RealBasis | BlochBasis | EigenBasis:
    """Construct one of :class:`RealBasis`, :class:`BlochBasis`, or :class:`EigenBasis` by name."""
    key = basis_name.strip().lower()
    if key == "real":
        return RealBasis(
            model,
            supercell,
            disorder_w=disorder_w,
            seed=seed,
            disorder_seed=disorder_seed,
            disorder=disorder,
            csr_dir=csr_dir,
        )
    if key == "bloch":
        return BlochBasis(
            model,
            supercell,
            disorder_w=disorder_w,
            seed=seed,
            disorder_seed=disorder_seed,
            disorder=disorder,
        )
    if key == "eigen":
        return EigenBasis(model, supercell, disorder_w=disorder_w, seed=seed)
    raise ValueError(
        f"make_basis: unknown basis '{basis_name}' (expected 'real', 'bloch', or 'eigen')"
    )
