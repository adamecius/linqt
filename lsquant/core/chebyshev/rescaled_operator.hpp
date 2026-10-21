/**
 * @file rescaled_operator.hpp
 * @brief LinearOperator affinely rescaling a Hamiltonian into `[-1, 1]`.
 *
 * Precomputes `H~ = window.Scale() * H + window.Shift() * I` via `SparseMatrix::Rescaled`
 * when constructed from a `SparseMatrix` (reproducing the exact IEEE-754 rounding of the
 * legacy in-place `Rescale` path), or evaluates `a * (scale * H + shift * I) * x + b * y`
 * on the fly when constructed from a general matrix-free `LinearOperator`.
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
 * `[-1, 1]` without mutating `H`. When constructed from a `SparseMatrix`, holds its own
 * compressed `SparseMatrix` built via `H.Rescaled(window.Scale(), window.Shift())` for
 * bit-for-bit agreement with legacy CSR rescaling. When constructed from a general
 * `LinearOperator`, wraps `H` without copying and applies `Scale * H + Shift * I` on the fly.
 *
 * @warning Throws `std::invalid_argument` if `window.upper <= window.lower`.
 * @par Pitfall
 *      When constructed from a `LinearOperator`, the referenced operator `H` must outlive
 *      this `RescaledOperator` instance.
 */
class RescaledOperator : public LinearOperator {
public:
  /**
   * @brief Construct the rescaled operator `H~` from a `SparseMatrix` and `window`.
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
   * @brief Construct the rescaled operator `H~` from a general `LinearOperator` and `window`.
   *
   * @param H      Physical matrix-free Hamiltonian (must outlive `*this`).
   * @param window Spectral window `[lower, upper]` with `upper > lower`.
   *
   * @warning Throws `std::invalid_argument` if `!(window.upper > window.lower)`.
   * @par Pitfall
   *      Does not copy `H`; keep `H` alive for the lifetime of this `RescaledOperator`.
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
  const LinearOperator* op_ = nullptr;
  Index dimension_ = 0;
  bool isHermitian_ = false;
};

namespace chebyshev {
using lsquant::RescaledOperator;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_RESCALED_OPERATOR_HPP_
