/**
 * @file rescaled_operator.hpp
 * @brief Dimensionless rescaled operators mapping physical spectra into `[-1, 1]`.
 *
 * Provides `RescaledOperator` (either owning a precomputed `a * H + b * I` CSR copy or
 * evaluating `a * H * x + b * x` matrix-free on any `LinearOperator`) and
 * `GeneralisedRescaledOperator` (owning `a * H' + b * S'` in CSR form and applying an
 * inverse overlap operator on every multiplication).
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_RESCALED_OPERATOR_HPP_
#define LSQUANT_CORE_CHEBYSHEV_RESCALED_OPERATOR_HPP_

#include <optional>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/chebyshev/spectral_window.hpp"

namespace lsquant {

/**
 * @brief Dimensionless rescaled operator `H~ = Scale * H + Shift * I` for Chebyshev recursion.
 *
 * Solves the problem of mapping a physical Hamiltonian `H` into the Chebyshev domain
 * `[-1, 1]` without mutating `H`. When constructed from a `SparseMatrix`, holds its own
 * compressed `SparseMatrix` built via `H.Rescaled(window.Scale(), window.Shift())` to
 * maintain bit-for-bit agreement with legacy in-place CSR rescaling (`O(dim + nnz)` setup
 * and `O(nnz)` per `Multiply`). When constructed from a general `LinearOperator` (such as
 * a two-sided Löwdin sandwich `S^{-1/2} H' S^{-1/2}`), evaluates
 * `y = a * (Scale * H * x + Shift * x) + b * y` matrix-free without materialising a CSR
 * copy.
 *
 * @warning Throws `std::invalid_argument` if `window.upper <= window.lower`.
 * @par Pitfall
 *      When wrapping a `LinearOperator` reference, the referenced operator must outlive this
 *      `RescaledOperator` instance.
 */
class RescaledOperator : public LinearOperator {
public:
  /**
   * @brief Construct the rescaled operator `H~` from a sparse matrix `H` and `window`.
   *
   * @param H      Physical sparse Hamiltonian (not modified).
   * @param window Spectral window `[lower, upper]` with `upper > lower`.
   *
   * @warning Throws `std::invalid_argument` if `!(window.upper > window.lower)`.
   * @par Pitfall
   *      `window` must enclose the entire eigenvalue spectrum of `H` or the Chebyshev
   *      recursion will diverge.
   */
  RescaledOperator(const SparseMatrix& H, const SpectralWindow& window);

  /**
   * @brief Construct a matrix-free rescaled operator `H~` over an arbitrary `LinearOperator`.
   *
   * Evaluates `y = a * (window.Scale() * H * x + window.Shift() * x) + b * y` without
   * building a CSR copy.
   *
   * @param H      Physical linear operator (must outlive `*this`).
   * @param window Spectral window `[lower, upper]` with `upper > lower`.
   *
   * @warning Throws `std::invalid_argument` if `!(window.upper > window.lower)` or if
   *          `H.Dimension() < 0`.
   * @par Pitfall
   *      Use the `SparseMatrix` overload when `H` is already a standalone `SparseMatrix` so
   *      that diagonal shifts use the legacy CSR rounding order.
   */
  RescaledOperator(const LinearOperator& H, const SpectralWindow& window);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Compute `y = a * H~ * x + b * y`.
   *
   * @param a Scalar multiplier for `H~ * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return cached Hermiticity of the rescaled operator.
   *
   * @return `true` if `H~` is Hermitian, `false` otherwise.
   */
  bool IsHermitian() const override;

  /**
   * @brief Return the spectral window used to build `H~`.
   *
   * @return Copy of the `SpectralWindow`.
   */
  const SpectralWindow& Window() const { return window_; }

private:
  SpectralWindow window_;
  std::optional<SparseMatrix> rescaled_;
  const LinearOperator* op_ = nullptr;
  bool isHermitian_ = false;
};

/**
 * @brief Generalised rescaled operator `H~ = S^{-1} (Scale * H' + Shift * S')`.
 *
 * Solves the problem of mapping the generalised pencil `H' c = E S' c` into the Chebyshev
 * interval `[-1, 1]` without mutating `H'`. Holds its own compressed CSR copy of the pencil
 * combination `Scale * H' + Shift * S'` (reproducing legacy `Rescale_nonOrth` arithmetic)
 * and applies `inverse` (`S^{-1}`) on every `Multiply` call. Construction takes
 * `O(dim + nnz(H') + nnz(S'))` time and space; each `Multiply` takes
 * `O(nnz(H' + S') + cost(S^{-1}))` time.
 *
 * @warning Always reports `IsHermitian() == false` because `S^{-1} H'` is self-adjoint
 *          with respect to `S'` rather than the Euclidean inner product.
 * @par Pitfall
 *      Never rescale as `S^{-1} (Scale * H' + Shift * I)`; shifting by `I` instead of `S'`
 *      adds `Shift * S^{-1}` rather than `Shift * I` to `Scale * S^{-1} H'`, corrupting
 *      every Chebyshev polynomial of order `m >= 1`.
 */
class GeneralisedRescaledOperator : public LinearOperator {
public:
  /**
   * @brief Construct `S^{-1} (window.Scale() * H + window.Shift() * S)`.
   *
   * @param H       Sparse Hamiltonian matrix `H'` (not modified).
   * @param S       Sparse overlap matrix `S'` (not modified; must match `H.Dimension()`).
   * @param inverse Linear operator applying `S^{-1}` (must match `H.Dimension()` and outlive `*this`).
   * @param window  Spectral window `[lower, upper]` with `upper > lower`.
   *
   * @warning Throws `std::invalid_argument` if dimensions differ or if
   *          `!(window.upper > window.lower)`.
   * @par Pitfall
   *      `inverse` is stored by const pointer and must remain valid for the lifetime of this
   *      instance.
   */
  GeneralisedRescaledOperator(const SparseMatrix& H, const SparseMatrix& S,
                              const LinearOperator& inverse, const SpectralWindow& window);

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Compute `y = a * S^{-1} * ((Scale * H' + Shift * S') * x) + b * y`.
   *
   * @param a Scalar multiplier for `H~ * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return `false` because `S^{-1} (Scale * H' + Shift * S')` is not Euclidean-Hermitian.
   *
   * @return Always `false`.
   */
  bool IsHermitian() const override { return false; }

  /**
   * @brief Return the spectral window used to build `H~`.
   *
   * @return Copy of the `SpectralWindow`.
   */
  const SpectralWindow& Window() const { return window_; }

private:
  SpectralWindow window_;
  SparseMatrix pencil_;
  const LinearOperator* inverse_ = nullptr;
};

namespace chebyshev {
using lsquant::GeneralisedRescaledOperator;
using lsquant::RescaledOperator;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_RESCALED_OPERATOR_HPP_
