"""Controlled non-orthogonal basis transformations for tight-binding CSR fixtures.

Provides ``make_local_transform`` to build a sparse, narrow-banded basis change
``R = I + alpha * K`` sharing the off-diagonal sparsity graph of a Hamiltonian ``H``,
and ``rotate_basis`` to map any dictionary of orthogonal CSR operators into the
transformed representation ``S' = R^+ R`` and ``O' = R^+ O R``. Also provides ASCII
``.CSR`` reader and writer routines matching ``lsquant::ReadCsrText`` and
``lsquant::WriteCsrText``.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp


def hermitian_fold(matrix: sp.csr_matrix) -> sp.csr_matrix:
    """Return ``0.5 * (M + M^+)`` in canonical sorted CSR layout.

    Reproduces ``lsquant::SparseMatrix::HermitianFold`` (and legacy
    ``SparseMatrixType::ConvertFromCSR``) so that operators read in Python match
    the C++ ``ReadCsrText`` representation down to IEEE-754 rounding.

    Parameters
    ----------
    matrix : scipy.sparse.csr_matrix
        Square complex sparse matrix.

    Returns
    -------
    scipy.sparse.csr_matrix
        Hermitian-folded CSR matrix with sorted column indices.
    """
    csr = matrix.tocsr().astype(np.complex128, copy=False)
    folded = (0.5 * (csr.conj().T + csr)).tocsr()
    folded.sum_duplicates()
    folded.sort_indices()
    return folded


def read_csr_text(path: str | Path, fold_hermitian: bool = True) -> sp.csr_matrix:
    """Read an LSQUANT ASCII ``.CSR`` operator file into a ``scipy.sparse.csr_matrix``.

    Parses ``dim nnz``, followed by ``nnz`` pairs of ``(re, im)`` floating-point values,
    ``nnz`` zero-based column indices, and ``dim + 1`` row pointer offsets.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the ``.CSR`` file.
    fold_hermitian : bool, default True
        When True, applies ``0.5 * (M^+ + M)`` to match ``lsquant::ReadCsrText``.

    Returns
    -------
    scipy.sparse.csr_matrix
        Complex CSR matrix of shape ``(dim, dim)``.
    """
    text = Path(path).read_text(encoding="utf-8").split()
    if len(text) < 2:
        raise ValueError(f"read_csr_text: empty or truncated header in {path}")
    dim = int(text[0])
    nnz = int(text[1])
    expected_tokens = 2 + 2 * nnz + nnz + (dim + 1)
    if len(text) < expected_tokens:
        raise ValueError(
            f"read_csr_text: expected {expected_tokens} tokens in {path}, got {len(text)}"
        )

    val_tokens = np.array(text[2 : 2 + 2 * nnz], dtype=np.float64)
    values = val_tokens[0::2] + 1j * val_tokens[1::2]

    col_offset = 2 + 2 * nnz
    columns = np.array(text[col_offset : col_offset + nnz], dtype=np.int64)

    rp_offset = col_offset + nnz
    row_pointers = np.array(text[rp_offset : rp_offset + dim + 1], dtype=np.int64)

    mat = sp.csr_matrix((values, columns, row_pointers), shape=(dim, dim), dtype=np.complex128)
    mat.sort_indices()
    if fold_hermitian:
        return hermitian_fold(mat)
    return mat


def write_csr_text(path: str | Path, matrix: sp.csr_matrix) -> None:
    """Write a ``scipy.sparse.csr_matrix`` to an LSQUANT ASCII ``.CSR`` operator file.

    Serialises ``dim nnz``, the ``nnz`` complex pairs ``re im`` at 17-digit precision,
    the ``nnz`` zero-based column indices, and the ``dim + 1`` CSR row pointers.

    Parameters
    ----------
    path : str or pathlib.Path
        Destination ``.CSR`` file path.
    matrix : scipy.sparse.csr_matrix
        Square sparse matrix to serialise.
    """
    csr = matrix.tocsr().astype(np.complex128, copy=False)
    csr.sum_duplicates()
    csr.sort_indices()
    if csr.shape[0] != csr.shape[1]:
        raise ValueError("write_csr_text: matrix must be square")

    dim = int(csr.shape[0])
    nnz = int(csr.nnz)
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Normalise -0.0 to +0.0 for deterministic ASCII output.
    re_vals = np.where(csr.data.real == 0.0, 0.0, csr.data.real)
    im_vals = np.where(csr.data.imag == 0.0, 0.0, csr.data.imag)

    val_strs = [f"{re:.17g} {im:.17g}" for re, im in zip(re_vals, im_vals)]
    col_strs = [str(int(c)) for c in csr.indices]
    rp_strs = [str(int(r)) for r in csr.indptr]

    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(f"{dim} {nnz}\n")
        fh.write(" ".join(val_strs) + (" \n" if val_strs else "\n"))
        fh.write(" ".join(col_strs) + (" \n" if col_strs else "\n"))
        fh.write(" ".join(rp_strs) + " \n")


def make_local_transform(
    H: sp.csr_matrix,
    alpha: float = 0.55,
    complex_entries: bool = True,
    seed: int = 12345,
) -> sp.csr_matrix:
    """Build a local basis transformation ``R = I + alpha * K`` on the sparsity graph of ``H``.

    Off-diagonal entries of ``K`` are drawn deterministically from
    ``numpy.random.default_rng(seed)`` (either real non-symmetric when
    ``complex_entries=False`` or complex non-Hermitian ``K_r + i K_i`` when
    ``complex_entries=True``) and normalised so that ``||K||_2 == 1`` and
    ``||alpha * K||_2 <= 0.55``. For ``N <= 4096``, verifies via dense ``eigvalsh``
    that ``S' = R^+ R`` is positive definite and ``cond(S')`` lies in ``[5.0, 10.0]``,
    raising ``ValueError`` outside that band.

    Parameters
    ----------
    H : scipy.sparse.csr_matrix
        Reference Hamiltonian defining the dimension ``N`` and off-diagonal sparsity graph.
    alpha : float, default 0.55
        Coupling amplitude multiplying ``K`` in ``R = I + alpha * K``.
    complex_entries : bool, default True
        If True, generate complex non-Hermitian off-diagonal entries; if False, generate
        real non-symmetric entries.
    seed : int, default 12345
        Seed for ``numpy.random.default_rng``.

    Returns
    -------
    scipy.sparse.csr_matrix
        Invertible transformation matrix ``R = I + alpha * K`` in CSR format.
    """
    csr = H.tocsr()
    if csr.shape[0] != csr.shape[1]:
        raise ValueError("make_local_transform: H must be square")
    dim = int(csr.shape[0])
    coo = csr.tocoo()
    mask = coo.row != coo.col
    rows = coo.row[mask]
    cols = coo.col[mask]
    nnz_off = int(rows.size)

    rng = np.random.default_rng(seed)
    if complex_entries:
        kr = rng.uniform(-1.0, 1.0, size=nnz_off)
        ki = rng.uniform(-1.0, 1.0, size=nnz_off)
        vals = (kr + 1j * ki).astype(np.complex128)
    else:
        kr = rng.uniform(-1.0, 1.0, size=nnz_off)
        vals = kr.astype(np.complex128)

    K = sp.csr_matrix((vals, (rows, cols)), shape=(dim, dim), dtype=np.complex128)
    K.sum_duplicates()
    K.sort_indices()

    if nnz_off > 0 and alpha != 0.0:
        sigma_max = float(np.linalg.norm(K.toarray(), ord=2))
        if sigma_max > 0.0:
            target_k_norm = min(1.0, 0.55 / abs(float(alpha)))
            K = K * (target_k_norm / sigma_max)

    ident = sp.eye(dim, format="csr", dtype=np.complex128)
    R = (ident + float(alpha) * K).tocsr()
    R.sum_duplicates()
    R.sort_indices()

    if dim <= 4096:
        s_mat = hermitian_fold((R.conj().T @ R).tocsr())
        s_eigs = la.eigvalsh(s_mat.toarray(), driver="evd")
        lmin = float(np.min(s_eigs))
        lmax = float(np.max(s_eigs))
        if lmin <= 0.0:
            raise ValueError(
                f"make_local_transform: S' is not positive definite (lambda_min={lmin})"
            )
        cond_s = lmax / lmin
        if not (5.0 <= cond_s <= 10.0):
            raise ValueError(
                f"make_local_transform: cond(S') = {cond_s:.6f} outside [5.0, 10.0]"
            )

    return R


def rotate_basis(
    ops: dict[str, sp.csr_matrix],
    R: sp.csr_matrix,
) -> dict[str, sp.csr_matrix]:
    """Return ``{'S': R.H @ R, **{k: R.H @ O @ R for k, O in ops.items()}}``.

    Applies the congruent transformation ``O' = R^+ O R`` to every operator in ``ops``
    and forms the overlap metric ``S' = R^+ R``.

    Parameters
    ----------
    ops : dict[str, scipy.sparse.csr_matrix]
        Dictionary mapping operator names (e.g. ``'HAM'``, ``'VX'``, ``'VY'``) to
        orthogonal-basis CSR matrices.
    R : scipy.sparse.csr_matrix
        Basis transformation matrix ``R`` of shape ``(N, N)``.

    Returns
    -------
    dict[str, scipy.sparse.csr_matrix]
        Dictionary containing ``'S'`` and every transformed operator in canonical CSR form.
    """
    r_csr = R.tocsr().astype(np.complex128, copy=False)
    r_dag = r_csr.conj().T.tocsr()

    s_mat = hermitian_fold((r_dag @ r_csr).tocsr())
    result: dict[str, sp.csr_matrix] = {"S": s_mat}
    for key, op in ops.items():
        op_csr = op.tocsr().astype(np.complex128, copy=False)
        rotated = (r_dag @ op_csr @ r_csr).tocsr()
        result[key] = hermitian_fold(rotated)
    return result
