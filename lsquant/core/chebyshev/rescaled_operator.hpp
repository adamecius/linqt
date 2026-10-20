/**
 * @file rescaled_operator.hpp
 * @brief LinearOperator holding an owned, affinely rescaled CSR copy of a Hamiltonian.
 *
 * Precomputes `H~ = window.Scale() * H + window.Shift() * I` via `SparseMatrix::Rescaled`
 * without mutating the caller's matrix `H`, reproducing the exact IEEE-754 rounding of the
 * legacy in-place `Rescale` path.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_RESCALED_OPERATOR_HPP_
#define LSQUANT_CORE_CHEBYSHEV_RESCALED_OPERATOR_HPP_

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/chebyshev/spectral_window.hpp"

namespace lsquant {

/**
 * @brief Dimensionless rescaled operator `H~ = Scale * H + Shift * I` for Chebyshev recursion.
 *
 * Solves the problem of mapping a physical Hamiltonian `H` into the Chebyshev domain
 * `[-1, 1]` without mutating `H` while maintaining bit-for-bit agreement with the legacy
 * in-place CSR rescaling. Holds its own compressed `SparseMatrix` built via
 * `H.Rescaled(window.Scale(), window.Shift())`. Construction takes `O(dim + nnz)` time and
 * space; each `Multiply` takes `O(nnz)` time.
 *
 * @warning Throws `std::invalid_argument` if `window.upper <= window.lower`.
 * @par Pitfall
 *      Evaluating `Scale * (H * x) + Shift * x` on the fly instead of multiplying by the
 *      rescaled CSR matrix changes floating-point rounding and breaks byte identity with
 *      legacy goldens.
 */
class RescaledOperator : public LinearOperator {
public:
  /**
   * @brief Construct the rescaled operator `H~` from `H` and `window`.
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
   * @brief Return the operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Compute `y = a * H~ * x + b * y` using the owned rescaled CSR matrix.
   *
   * @param a Scalar multiplier for `H~ * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Scalar multiplier for `y`.
   * @param y Input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return cached Hermiticity of the rescaled matrix.
   *
   * @return `true` if `H~` is Hermitian within `1e-12`, `false` otherwise.
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
  SparseMatrix rescaled_;
  bool isHermitian_;
};

namespace chebyshev {
using lsquant::RescaledOperator;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_RESCALED_OPERATOR_HPP_
