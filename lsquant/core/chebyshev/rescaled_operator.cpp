/**
 * @file rescaled_operator.cpp
 * @brief Implementation of RescaledOperator.
 */
#include "core/chebyshev/rescaled_operator.hpp"

#include <cmath>
#include <stdexcept>

namespace lsquant {

namespace {

void ValidateWindow(const SpectralWindow& window) {
  if (!(window.upper > window.lower) || !std::isfinite(window.lower) ||
      !std::isfinite(window.upper)) {
    throw std::invalid_argument(
        "RescaledOperator: window.upper must be greater than window.lower");
  }
}

}  // namespace

RescaledOperator::RescaledOperator(const SparseMatrix& H, const SpectralWindow& window)
    : window_(window) {
  ValidateWindow(window);
  rescaled_ = H.Rescaled(Scalar(window.Scale(), 0.0), Scalar(window.Shift(), 0.0));
  dimension_ = rescaled_.Dimension();
  isHermitian_ = rescaled_.IsHermitian();
}

RescaledOperator::RescaledOperator(const LinearOperator& H, const SpectralWindow& window)
    : window_(window),
      op_(&H),
      dimension_(H.Dimension()),
      isHermitian_(H.IsHermitian()) {
  ValidateWindow(window);
}

Index RescaledOperator::Dimension() const { return dimension_; }

void RescaledOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  if (op_ == nullptr) {
    rescaled_.Multiply(a, x, b, y);
    return;
  }

  const Index dim = dimension_;
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("RescaledOperator::Multiply: vector dimension mismatch");
  }

  const Scalar aScale = a * Scalar(window_.Scale(), 0.0);
  const Scalar aShift = a * Scalar(window_.Shift(), 0.0);
  if (&x == &y) {
    Vector xCopy = x;
    op_->Multiply(aScale, xCopy, b, y);
    if (aShift != Scalar(0.0, 0.0)) {
      const Scalar* xPtr = xCopy.data();
      Scalar* yPtr = y.data();
      #pragma omp parallel for schedule(static)
      for (Index i = 0; i < dim; ++i) {
        yPtr[i] += aShift * xPtr[i];
      }
    }
    return;
  }

  op_->Multiply(aScale, x, b, y);
  if (aShift != Scalar(0.0, 0.0)) {
    const Scalar* xPtr = x.data();
    Scalar* yPtr = y.data();
    #pragma omp parallel for schedule(static)
    for (Index i = 0; i < dim; ++i) {
      yPtr[i] += aShift * xPtr[i];
    }
  }
}

bool RescaledOperator::IsHermitian() const { return isHermitian_; }

}  // namespace lsquant
