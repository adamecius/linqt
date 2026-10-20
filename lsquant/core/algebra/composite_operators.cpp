/**
 * @file composite_operators.cpp
 * @brief Implementation of algebraic composite operators.
 */
#include "core/algebra/composite_operators.hpp"

#include <stdexcept>

namespace lsquant {

IdentityOperator::IdentityOperator(Index dimension) : dimension_(dimension) {
  if (dimension < 0) {
    throw std::invalid_argument("IdentityOperator: negative dimension");
  }
}

Index IdentityOperator::Dimension() const { return dimension_; }

void IdentityOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  if (static_cast<Index>(x.size()) != dimension_ ||
      static_cast<Index>(y.size()) != dimension_) {
    throw std::invalid_argument("IdentityOperator::Multiply: vector dimension mismatch");
  }
  for (Index i = 0; i < dimension_; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    y[idx] = a * x[idx] + b * y[idx];
  }
}

ScaledOperator::ScaledOperator(Scalar alpha, const LinearOperator& op)
    : alpha_(alpha), op_(&op) {}

Index ScaledOperator::Dimension() const { return op_->Dimension(); }

void ScaledOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  op_->Multiply(a * alpha_, x, b, y);
}

bool ScaledOperator::IsHermitian() const {
  return op_->IsHermitian() && (alpha_.imag() == 0.0);
}

OperatorSum::OperatorSum(const LinearOperator& first, const LinearOperator& second)
    : first_(&first), second_(&second) {
  if (first.Dimension() != second.Dimension()) {
    throw std::invalid_argument("OperatorSum: operator dimensions must match");
  }
}

Index OperatorSum::Dimension() const { return first_->Dimension(); }

void OperatorSum::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("OperatorSum::Multiply: vector dimension mismatch");
  }
  if (&x == &y) {
    Vector xCopy = x;
    first_->Multiply(a, xCopy, b, y);
    second_->Multiply(a, xCopy, Scalar(1.0, 0.0), y);
    return;
  }
  first_->Multiply(a, x, b, y);
  second_->Multiply(a, x, Scalar(1.0, 0.0), y);
}

bool OperatorSum::IsHermitian() const {
  return first_->IsHermitian() && second_->IsHermitian();
}

OperatorProduct::OperatorProduct(const LinearOperator& first, const LinearOperator& second)
    : first_(&first), second_(&second) {
  if (first.Dimension() != second.Dimension()) {
    throw std::invalid_argument("OperatorProduct: operator dimensions must match");
  }
}

Index OperatorProduct::Dimension() const { return first_->Dimension(); }

void OperatorProduct::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("OperatorProduct::Multiply: vector dimension mismatch");
  }
  Vector temp(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  second_->Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), temp);
  first_->Multiply(a, temp, b, y);
}

bool OperatorProduct::IsHermitian() const {
  if (first_->IsIdentity()) {
    return second_->IsHermitian();
  }
  if (second_->IsIdentity()) {
    return first_->IsHermitian();
  }
  return false;
}

AdjointOperator::AdjointOperator(const LinearOperator& op) : isHermitian_(op.IsHermitian()) {
  if (const auto* sparseOp = dynamic_cast<const SparseMatrixOperator*>(&op)) {
    adjointMatrix_ = sparseOp->Matrix().Adjoint();
  } else if (op.IsHermitian()) {
    forwardedOp_ = &op;
  } else {
    throw std::invalid_argument(
        "AdjointOperator: non-Hermitian operator must be a SparseMatrixOperator");
  }
}

AdjointOperator::AdjointOperator(const SparseMatrix& matrix)
    : adjointMatrix_(matrix.Adjoint()), isHermitian_(matrix.IsHermitian()) {}

Index AdjointOperator::Dimension() const {
  if (adjointMatrix_.has_value()) {
    return adjointMatrix_->Dimension();
  }
  return forwardedOp_->Dimension();
}

void AdjointOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  if (adjointMatrix_.has_value()) {
    adjointMatrix_->Multiply(a, x, b, y);
  } else {
    forwardedOp_->Multiply(a, x, b, y);
  }
}

bool AdjointOperator::IsHermitian() const { return isHermitian_; }

LowdinOperator::LowdinOperator(const LinearOperator& inverseSqrt, const LinearOperator& O)
    : inverseSqrt_(&inverseSqrt), op_(&O) {
  if (inverseSqrt.Dimension() != O.Dimension()) {
    throw std::invalid_argument("LowdinOperator: operator dimensions must match");
  }
}

Index LowdinOperator::Dimension() const { return op_->Dimension(); }

void LowdinOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("LowdinOperator::Multiply: vector dimension mismatch");
  }
  Vector t1(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector t2(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  inverseSqrt_->Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), t1);
  op_->Multiply(Scalar(1.0, 0.0), t1, Scalar(0.0, 0.0), t2);
  inverseSqrt_->Multiply(a, t2, b, y);
}

bool LowdinOperator::IsHermitian() const {
  return inverseSqrt_->IsHermitian() && op_->IsHermitian();
}

MetricInverseProduct::MetricInverseProduct(const LinearOperator& inverse,
                                           const LinearOperator& O)
    : inverse_(&inverse), op_(&O) {
  if (inverse.Dimension() != O.Dimension()) {
    throw std::invalid_argument("MetricInverseProduct: operator dimensions must match");
  }
}

Index MetricInverseProduct::Dimension() const { return op_->Dimension(); }

void MetricInverseProduct::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("MetricInverseProduct::Multiply: vector dimension mismatch");
  }
  Vector temp(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  op_->Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), temp);
  inverse_->Multiply(a, temp, b, y);
}

}  // namespace lsquant
