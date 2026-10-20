"""Type stubs for the pybind11 module ``lsquant._core``."""

from __future__ import annotations

from typing import Sequence, Union, overload
import numpy as np
import numpy.typing as npt

class SparseMatrix:
    def __init__(self) -> None: ...
    @staticmethod
    def from_csr(
        dim: int,
        row_pointers: npt.NDArray[np.int64],
        columns: npt.NDArray[np.int64],
        values: npt.NDArray[np.complex128],
    ) -> SparseMatrix: ...
    @property
    def dimension(self) -> int: ...
    @property
    def non_zeros(self) -> int: ...
    def hermitian_fold(self) -> SparseMatrix: ...
    def is_hermitian(self, tol: float = 1e-12) -> bool: ...

class LinearOperator:
    @property
    def dimension(self) -> int: ...
    @property
    def is_identity(self) -> bool: ...
    @property
    def is_hermitian(self) -> bool: ...

class SparseMatrixOperator(LinearOperator):
    def __init__(self, matrix: SparseMatrix) -> None: ...

class IdentityOperator(LinearOperator):
    def __init__(self, dim: int) -> None: ...

class LowdinOperator(LinearOperator):
    def __init__(self, inverse_sqrt: LinearOperator, operator: LinearOperator) -> None: ...

class MetricInverseProduct(LinearOperator):
    def __init__(self, inverse: LinearOperator, operator: LinearOperator) -> None: ...

class MetricPolicy:
    @property
    def is_orthogonal(self) -> bool: ...
    @property
    def overlap(self) -> LinearOperator: ...
    @property
    def inverse_overlap(self) -> LinearOperator: ...
    @property
    def inverse_sqrt_overlap(self) -> LinearOperator: ...

class OrthogonalMetric(MetricPolicy):
    def __init__(self, dim: int = 0) -> None: ...

class NonOrthogonalMetric(MetricPolicy):
    def __init__(
        self, s: SparseMatrix, inverse: LinearOperator, inverse_sqrt: LinearOperator
    ) -> None: ...

class ExactInverseOverlap(LinearOperator):
    def __init__(self, s: SparseMatrix, max_dimension: int = 4096) -> None: ...

class ExactInverseSqrtOverlap(LinearOperator):
    def __init__(self, s: SparseMatrix, max_dimension: int = 4096) -> None: ...

class IschoInverseOverlap(LinearOperator):
    def __init__(
        self, s: SparseMatrix, tolerance: float = 1e-13, max_iterations: int = 200
    ) -> None: ...
    @property
    def last_iterations(self) -> int: ...
    @property
    def last_error(self) -> float: ...

class ChebyshevInverseOverlap(LinearOperator):
    def __init__(
        self,
        s: SparseMatrix,
        bounds_of_s: SpectralBounds,
        order: int,
        coeff_tol: float = 1e-12,
    ) -> None: ...
    @property
    def effective_order(self) -> int: ...

class ChebyshevInverseSqrtOverlap(LinearOperator):
    def __init__(
        self,
        s: SparseMatrix,
        bounds_of_s: SpectralBounds,
        order: int,
        coeff_tol: float = 1e-12,
    ) -> None: ...
    @property
    def effective_order(self) -> int: ...

class StateGenerator:
    @property
    def dimension(self) -> int: ...
    @property
    def number_of_states(self) -> int: ...
    @property
    def label(self) -> str: ...
    def reset(self) -> None: ...

class RandomPhaseStates(StateGenerator):
    def __init__(self, dim: int, num_states: int, seed: int = 12345) -> None: ...

class LocalStates(StateGenerator):
    def __init__(self, dim: int, sites: Union[npt.NDArray[np.int64], Sequence[int]]) -> None: ...

class UserStates(StateGenerator):
    def __init__(
        self, states: Sequence[Union[npt.NDArray[np.complex128], Sequence[complex]]]
    ) -> None: ...

class SpectralBounds:
    lower: float
    upper: float
    def __init__(self, lower: float = 0.0, upper: float = 0.0) -> None: ...

def GershgorinBounds(matrix: SparseMatrix, pad: float = 0.10) -> SpectralBounds: ...
def LanczosBounds(
    operator: LinearOperator, iterations: int = 40, seed: int = 12345
) -> SpectralBounds: ...

class SpectralWindow:
    lower: float
    upper: float
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, lower: float, upper: float) -> None: ...
    @property
    def scale(self) -> float: ...
    @property
    def shift(self) -> float: ...
    @property
    def half_width(self) -> float: ...
    @property
    def center(self) -> float: ...

class RescaledOperator(LinearOperator):
    @overload
    def __init__(self, matrix: SparseMatrix, window: SpectralWindow) -> None: ...
    @overload
    def __init__(self, operator: LinearOperator, window: SpectralWindow) -> None: ...
    @property
    def window(self) -> SpectralWindow: ...

class GeneralisedRescaledOperator(LinearOperator):
    def __init__(
        self,
        h: SparseMatrix,
        s: SparseMatrix,
        inverse: LinearOperator,
        window: SpectralWindow,
    ) -> None: ...
    @property
    def window(self) -> SpectralWindow: ...

class ApplyOperator:
    operator_index: int
    def __init__(self, operator_index: int) -> None: ...

class Expand:
    num_moments: int
    def __init__(self, num_moments: int) -> None: ...

class ApplySeries:
    coeffs: list[complex]
    def __init__(
        self, coeffs: Union[npt.NDArray[np.complex128], Sequence[complex]]
    ) -> None: ...

class MirrorToBra:
    def __init__(self) -> None: ...

class Trace:
    weights: list[complex]
    def __init__(
        self, weights: Union[npt.NDArray[np.complex128], Sequence[complex]]
    ) -> None: ...

Instruction = Union[ApplyOperator, Expand, ApplySeries, MirrorToBra, Trace]

class Plan:
    code: list[Instruction]
    word: str
    def __init__(self, code: Sequence[Instruction] = ..., word: str = "") -> None: ...

class ChebyshevMoments:
    window: SpectralWindow
    system_size: int
    number_of_states: int
    state_label: str
    word: str
    kpm_weighted: bool
    symmetrised: bool
    def __init__(self, dims: Sequence[int] = ...) -> None: ...
    @property
    def rank(self) -> int: ...
    @property
    def dims(self) -> list[int]: ...
    @property
    def data(self) -> npt.NDArray[np.complex128]: ...

class MomentEngine:
    def __init__(
        self,
        rescaled_h: LinearOperator,
        operators: Sequence[LinearOperator],
        states: StateGenerator,
        batch_size: int,
    ) -> None: ...
    def run(self, plan: Plan) -> ChebyshevMoments: ...

def ApplyKpmWeights(
    moments: ChebyshevMoments, num_states: int, symmetrise: bool
) -> None: ...
def ReadCsrText(path: str) -> SparseMatrix: ...
def ReadChebmom(path: str) -> ChebyshevMoments: ...
def WriteChebmom(path: str, moments: ChebyshevMoments) -> None: ...
def WriteLegacyChebmom1D(path: str, moments: ChebyshevMoments) -> None: ...
def WriteLegacyChebmom2D(path: str, moments: ChebyshevMoments) -> None: ...

@overload
def ContractRank1(
    moments: npt.NDArray[np.float64], table: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]: ...
@overload
def ContractRank1(
    moments: ChebyshevMoments, table: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]: ...

@overload
def ContractRank2Imag(
    moments: npt.NDArray[np.complex128],
    table_a: npt.NDArray[np.float64],
    table_b: npt.NDArray[np.complex128],
) -> npt.NDArray[np.float64]: ...
@overload
def ContractRank2Imag(
    moments: ChebyshevMoments,
    table_a: npt.NDArray[np.float64],
    table_b: npt.NDArray[np.complex128],
) -> npt.NDArray[np.float64]: ...

@overload
def ContractRank2Real(
    moments: npt.NDArray[np.complex128],
    table_a: npt.NDArray[np.complex128],
    table_b: npt.NDArray[np.complex128],
) -> npt.NDArray[np.float64]: ...
@overload
def ContractRank2Real(
    moments: ChebyshevMoments,
    table_a: npt.NDArray[np.complex128],
    table_b: npt.NDArray[np.complex128],
) -> npt.NDArray[np.float64]: ...

@overload
def ContractRank2RealReal(
    moments: npt.NDArray[np.complex128],
    table_a: npt.NDArray[np.float64],
    table_b: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]: ...
@overload
def ContractRank2RealReal(
    moments: ChebyshevMoments,
    table_a: npt.NDArray[np.float64],
    table_b: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]: ...

def WriteCurveText(
    path: str, x: npt.NDArray[np.float64], y: npt.NDArray[np.float64]
) -> None: ...

class compat:
    @staticmethod
    def legacy_green_tables(
        x_grid: Union[Sequence[float], npt.NDArray[np.float64]],
        num_moments_0: int,
        num_moments_1: int | None = None,
    ) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]: ...
