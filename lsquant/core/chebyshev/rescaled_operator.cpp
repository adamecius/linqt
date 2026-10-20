/**
 * @file rescaled_operator.cpp
 * @brief Implementation of RescaledOperator.
 */
#include "core/chebyshev/rescaled_operator.hpp"

#include <cmath>
#include <stdexcept>

namespace lsquant {

RescaledOperator::RescaledOperator(const SparseMatrix& H, const SpectralWindow& window)
    : window_(window) {
  if (!(window.upper > window.lower) || !std::isfinite(window.lower) ||
      !std::isfinite(window.upper)) {
    throw std::invalid_argument("RescaledOperator: window.upper must be greater than window.lower");
  }
  rescaled_ = H.Rescaled(Scalar(window.Scale(), 0.0), Scalar(window.Shift(), 0.0));
  isHermitian_ = rescaled_.IsHermitian();
}

Index RescaledOperator::Dimension() const { return rescaled_.Dimension(); }

void RescaledOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  rescaled_.Multiply(a, x, b, y);
}

bool RescaledOperator::IsHermitian() const { return isHermitian_; }

}  // namespace lsquant
