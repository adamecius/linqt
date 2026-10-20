/**
 * @file rescaled_operator.cpp
 * @brief Implementation of RescaledOperator and GeneralisedRescaledOperator.
 */
#include "core/chebyshev/rescaled_operator.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace lsquant {

namespace {

void ValidateWindow(const SpectralWindow& window, const char* caller) {
  if (!(window.upper > window.lower) || !std::isfinite(window.lower) ||
      !std::isfinite(window.upper)) {
    throw std::invalid_argument(std::string(caller) +
                                ": window.upper must be greater than window.lower");
  }
}

}  // namespace

RescaledOperator::RescaledOperator(const SparseMatrix& H, const SpectralWindow& window)
    : window_(window) {
  ValidateWindow(window, "RescaledOperator");
  rescaled_ = H.Rescaled(Scalar(window.Scale(), 0.0), Scalar(window.Shift(), 0.0));
  isHermitian_ = rescaled_->IsHermitian();
}

RescaledOperator::RescaledOperator(const LinearOperator& H, const SpectralWindow& window)
    : window_(window), op_(&H), isHermitian_(H.IsHermitian()) {
  ValidateWindow(window, "RescaledOperator");
  if (H.Dimension() < 0) {
    throw std::invalid_argument("RescaledOperator: negative operator dimension");
  }
}

Index RescaledOperator::Dimension() const {
  if (rescaled_.has_value()) {
    return rescaled_->Dimension();
  }
  return op_->Dimension();
}

void RescaledOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  if (rescaled_.has_value()) {
    rescaled_->Multiply(a, x, b, y);
    return;
  }
  const Index dim = op_->Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("RescaledOperator::Multiply: vector dimension mismatch");
  }
  Vector hx(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  op_->Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), hx);
  const Scalar scale(window_.Scale(), 0.0);
  const Scalar shift(window_.Shift(), 0.0);
  for (Index i = 0; i < dim; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    y[idx] = a * (scale * hx[idx] + shift * x[idx]) + b * y[idx];
  }
}

bool RescaledOperator::IsHermitian() const { return isHermitian_; }

GeneralisedRescaledOperator::GeneralisedRescaledOperator(const SparseMatrix& H,
                                                         const SparseMatrix& S,
                                                         const LinearOperator& inverse,
                                                         const SpectralWindow& window)
    : window_(window), inverse_(&inverse) {
  ValidateWindow(window, "GeneralisedRescaledOperator");
  const Index dim = H.Dimension();
  if (S.Dimension() != dim || inverse.Dimension() != dim) {
    throw std::invalid_argument("GeneralisedRescaledOperator: operator dimensions must match");
  }
  Eigen::SparseMatrix<Scalar, Eigen::RowMajor, Index> combined =
      Scalar(window.Scale(), 0.0) * H.Eigen() + Scalar(window.Shift(), 0.0) * S.Eigen();
  if (!combined.isCompressed()) {
    combined.makeCompressed();
  }
  const Index nnz = combined.nonZeros();
  const Index* const rpPtr = combined.outerIndexPtr();
  const Index* const clPtr = combined.innerIndexPtr();
  const Scalar* const vlPtr = combined.valuePtr();
  std::vector<Index> rowPointers(rpPtr, rpPtr + dim + 1);
  std::vector<Index> columns(clPtr, clPtr + nnz);
  std::vector<Scalar> values(vlPtr, vlPtr + nnz);
  pencil_ = SparseMatrix::FromCsr(dim, std::move(rowPointers), std::move(columns),
                                  std::move(values));
}

Index GeneralisedRescaledOperator::Dimension() const { return pencil_.Dimension(); }

void GeneralisedRescaledOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = pencil_.Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument(
        "GeneralisedRescaledOperator::Multiply: vector dimension mismatch");
  }
  Vector temp(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  pencil_.Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), temp);
  inverse_->Multiply(a, temp, b, y);
}

}  // namespace lsquant
