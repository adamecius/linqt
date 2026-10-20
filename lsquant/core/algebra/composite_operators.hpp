/**
 * @file composite_operators.hpp
 * @brief Algebraic composite operators implementing LinearOperator.
 *
 * Provides matrix-free operator combinators (`IdentityOperator`, `ScaledOperator`,
 * `OperatorSum`, `OperatorProduct`, `AdjointOperator`, `LowdinOperator`,
 * `MetricInverseProduct`) so that linear combinations, products, symmetric Löwdin
 * sandwiches, and one-sided metric-inverse products can be passed into the Chebyshev
 * engine without materialising intermediate matrices.
 */
#ifndef LSQUANT_CORE_ALGEBRA_COMPOSITE_OPERATORS_HPP_
#define LSQUANT_CORE_ALGEBRA_COMPOSITE_OPERATORS_HPP_

#include <optional>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Matrix-free identity operator `I` on a complex vector space of dimension `N`.
 *
 * Implements `y = a * x + b * y` in `O(N)` time and reports `IsIdentity() == true` so
 * that callers can bypass redundant identity applications in `O(1)`.
 *
 * @warning Throws `std::invalid_argument` if constructed with a negative dimension or
 *          applied to vectors of mismatched size.
 * @par Pitfall
 *      Even when `IsIdentity()` is true, calling `Multiply(a, x, b, y)` still applies
 *      the scalar weights `a` and `b`.
 */
class IdentityOperator : public LinearOperator {
public:
  /**
   * @brief Construct an identity operator of dimension `dimension`.
   *
   * @param dimension Vector-space dimension `N >= 0`.
   *
   * @warning Throws `std::invalid_argument` if `dimension < 0`.
   * @par Pitfall
   *      Ensure `dimension` matches the Hamiltonian dimension when mixing operators.
   */
  explicit IdentityOperator(Index dimension);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * x + b * y`.
   *
   * @param a Scalar multiplier for `x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` to indicate structural identity.
   *
   * @return Always `true`.
   */
  bool IsIdentity() const override { return true; }

  /**
   * @brief Return `true` because the identity map is self-adjoint.
   *
   * @return Always `true`.
   */
  bool IsHermitian() const override { return true; }

private:
  Index dimension_;
};

/**
 * @brief Matrix-free scalar multiple `alpha * O` of an existing `LinearOperator`.
 *
 * Forwards `Multiply(a, x, b, y)` to `op.Multiply(a * alpha, x, b, y)` with zero extra
 * vector allocations and `O(1)` overhead beyond the wrapped operator.
 *
 * @warning Does not own the wrapped `LinearOperator`; `op` must outlive this instance.
 * @par Pitfall
 *      A Hermitian operator scaled by a complex number with non-zero imaginary part is
 *      no longer Hermitian.
 */
class ScaledOperator : public LinearOperator {
public:
  /**
   * @brief Construct the scaled operator `alpha * op`.
   *
   * @param alpha Complex scalar factor.
   * @param op    Referenced operator (must outlive `*this`).
   */
  ScaledOperator(Scalar alpha, const LinearOperator& op);

  /**
   * @brief Return the dimension of the wrapped operator.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = (a * alpha) * O * x + b * y`.
   *
   * @param a Scalar multiplier for `(alpha * O) * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` iff the wrapped operator is Hermitian and `alpha` is real.
   *
   * @return Hermiticity of `alpha * O`.
   */
  bool IsHermitian() const override;

private:
  Scalar alpha_;
  const LinearOperator* op_;
};

/**
 * @brief Matrix-free sum `O1 + O2` of two linear operators of equal dimension.
 *
 * Evaluates `y = a * (O1 + O2) * x + b * y` as `y = a * O1 * x + b * y` followed by
 * `y = a * O2 * x + 1 * y`. Complexity is `O(cost(O1) + cost(O2))` with no temporary
 * vector allocation unless `x` and `y` alias.
 *
 * @warning Both referenced operators must outlive this `OperatorSum` instance.
 * @par Pitfall
 *      Throws `std::invalid_argument` at construction if `first.Dimension() != second.Dimension()`.
 */
class OperatorSum : public LinearOperator {
public:
  /**
   * @brief Construct the sum `first + second`.
   *
   * @param first  First summand operator.
   * @param second Second summand operator (must match `first.Dimension()`).
   *
   * @warning Throws `std::invalid_argument` if dimensions differ.
   * @par Pitfall
   *      Neither summand is copied; keep both alive for the lifetime of this object.
   */
  OperatorSum(const LinearOperator& first, const LinearOperator& second);

  /**
   * @brief Return the common operator dimension.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * (O1 + O2) * x + b * y`.
   *
   * @param a Scalar multiplier for `(O1 + O2) * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` if both summands report `IsHermitian() == true`.
   *
   * @return Sufficient condition for Hermiticity of `O1 + O2`.
   */
  bool IsHermitian() const override;

private:
  const LinearOperator* first_;
  const LinearOperator* second_;
};

/**
 * @brief Matrix-free composition `O1 * O2` of two linear operators (`O2` applied first).
 *
 * Evaluates `temp = O2 * x` into a temporary vector of length `Dimension()`, then
 * `y = a * O1 * temp + b * y`. Complexity is `O(cost(O1) + cost(O2) + N)`.
 *
 * @warning Both referenced operators must outlive this `OperatorProduct` instance.
 * @par Pitfall
 *      In operator notation `O1 * O2`, `second` (`O2`) acts on the ket `x` first, and
 *      `first` (`O1`) acts on the intermediate result.
 */
class OperatorProduct : public LinearOperator {
public:
  /**
   * @brief Construct the composition `first * second` (`second` applied first to `x`).
   *
   * @param first  Outer operator `O1`.
   * @param second Inner operator `O2` applied first to the input ket.
   *
   * @warning Throws `std::invalid_argument` if `first.Dimension() != second.Dimension()`.
   * @par Pitfall
   *      A product of two Hermitian operators is generally non-Hermitian unless they commute.
   */
  OperatorProduct(const LinearOperator& first, const LinearOperator& second);

  /**
   * @brief Return the common operator dimension.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * O1 * (O2 * x) + b * y`.
   *
   * @param a Scalar multiplier for `(O1 * O2) * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Report whether the product is guaranteed to be Hermitian.
   *
   * @return `true` if one factor is the identity and the other is Hermitian; `false` otherwise.
   */
  bool IsHermitian() const override;

private:
  const LinearOperator* first_;
  const LinearOperator* second_;
};

/**
 * @brief Conjugate transpose operator `O^+`.
 *
 * When constructed from a `SparseMatrix` or a `SparseMatrixOperator`, stores the explicit
 * adjoint CSR matrix (`SparseMatrix::Adjoint()`) so `Multiply` uses the row-parallel CSR
 * kernel even when `O` is non-Hermitian. When constructed from a general Hermitian
 * `LinearOperator`, forwards directly to `op.Multiply`.
 *
 * @warning Throws `std::invalid_argument` if constructed from a non-Hermitian matrix-free
 *          `LinearOperator` that is not backed by a `SparseMatrixOperator`.
 * @par Pitfall
 *      Materialising the adjoint CSR of a `SparseMatrixOperator` takes `O(dim + nnz)`
 *      memory and construction time; construct `AdjointOperator` outside inner loops.
 */
class AdjointOperator : public LinearOperator {
public:
  /**
   * @brief Construct the adjoint `op^+` of a `LinearOperator`.
   *
   * If `op` is a `SparseMatrixOperator`, builds and stores its adjoint CSR matrix.
   * Otherwise requires `op.IsHermitian() == true` and forwards to `op`.
   *
   * @param op Referenced operator.
   *
   * @warning Throws `std::invalid_argument` if `op` is neither a `SparseMatrixOperator`
   *          nor Hermitian.
   * @par Pitfall
   *      If `op` is forwarded without an owned CSR copy, `op` must outlive `*this`.
   */
  explicit AdjointOperator(const LinearOperator& op);

  /**
   * @brief Construct the adjoint `matrix^+` directly from a `SparseMatrix`.
   *
   * @param matrix Source sparse matrix whose conjugate transpose is stored.
   */
  explicit AdjointOperator(const SparseMatrix& matrix);

  /**
   * @brief Return the operator dimension.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * O^+ * x + b * y`.
   *
   * @param a Scalar multiplier for `O^+ * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Report whether `O^+` is Hermitian.
   *
   * @return `true` if `O` is Hermitian, `false` otherwise.
   */
  bool IsHermitian() const override;

private:
  std::optional<SparseMatrix> adjointMatrix_;
  const LinearOperator* forwardedOp_ = nullptr;
  bool isHermitian_ = false;
};

/**
 * @brief Matrix-free symmetric Löwdin sandwich `S^{-1/2} O S^{-1/2}`.
 *
 * Solves the problem of transforming a non-orthogonal basis representation `O'` of a
 * physical operator into the orthonormal Löwdin representation without materialising
 * dense `N x N` matrices. Evaluates `t1 = S^{-1/2} x`, `t2 = O t1`, and
 * `y = a * S^{-1/2} t2 + b * y` using two temporary vectors of length `Dimension()` per
 * `Multiply` call. Complexity is `O(2 * cost(S^{-1/2}) + cost(O) + N)`.
 *
 * @warning Neither `inverseSqrt` nor `O` is owned; both must outlive this instance.
 * @par Pitfall
 *      Both factors must have identical `Dimension()` or the constructor throws
 *      `std::invalid_argument`.
 */
class LowdinOperator : public LinearOperator {
public:
  /**
   * @brief Construct the symmetric sandwich `inverseSqrt * O * inverseSqrt`.
   *
   * @param inverseSqrt Inverse square-root overlap applier `S^{-1/2}`.
   * @param O           Inner operator `O'` in the non-orthogonal representation.
   *
   * @warning Throws `std::invalid_argument` if `inverseSqrt.Dimension() != O.Dimension()`.
   * @par Pitfall
   *      For `IsHermitian()` to return `true`, both `inverseSqrt` and `O` must report
   *      `IsHermitian() == true`.
   */
  LowdinOperator(const LinearOperator& inverseSqrt, const LinearOperator& O);

  /**
   * @brief Return the common operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * (S^{-1/2} * O * S^{-1/2}) * x + b * y`.
   *
   * @param a Scalar multiplier for `(S^{-1/2} O S^{-1/2}) x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `true` iff both `inverseSqrt` and `O` are Hermitian.
   *
   * @return Hermiticity of `S^{-1/2} O S^{-1/2}`.
   */
  bool IsHermitian() const override;

private:
  const LinearOperator* inverseSqrt_;
  const LinearOperator* op_;
};

/**
 * @brief Matrix-free one-sided product `S^{-1} O` for the generalised non-orthogonal route.
 *
 * Solves the problem of representing the effective operator `S^{-1} O'` in the
 * non-orthogonal basis where `O'` maps coordinates to covariant dual components and
 * `S^{-1}` raises the index back to contravariant coordinates. Evaluates `t = O x` into a
 * single temporary vector and then `y = a * S^{-1} t + b * y`. Complexity is
 * `O(cost(S^{-1}) + cost(O) + N)`.
 *
 * @warning Reports `IsHermitian() == false` unconditionally because `S^{-1} O` is
 *          `S`-self-adjoint rather than Euclidean-Hermitian whenever `S != I`.
 * @par Pitfall
 *      Never assume `(S^{-1} O)^+ = S^{-1} O`; the Euclidean adjoint is `O^+ S^{-1}`,
 *      which is why the generalised Kubo bra path requires explicit `S^{-1}` and `S`
 *      insertions.
 */
class MetricInverseProduct : public LinearOperator {
public:
  /**
   * @brief Construct the one-sided product `inverse * O`.
   *
   * @param inverse Inverse overlap applier `S^{-1}`.
   * @param O       Right factor `O'`.
   *
   * @warning Throws `std::invalid_argument` if `inverse.Dimension() != O.Dimension()`.
   * @par Pitfall
   *      Neither referenced operator is copied; both must outlive this object.
   */
  MetricInverseProduct(const LinearOperator& inverse, const LinearOperator& O);

  /**
   * @brief Return the common operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Evaluate `y = a * S^{-1} * (O * x) + b * y`.
   *
   * @param a Scalar multiplier for `(S^{-1} O) x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `false` because `S^{-1} O` is not Euclidean-Hermitian in general.
   *
   * @return Always `false`.
   */
  bool IsHermitian() const override { return false; }

private:
  const LinearOperator* inverse_;
  const LinearOperator* op_;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_COMPOSITE_OPERATORS_HPP_
