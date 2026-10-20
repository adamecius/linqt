/**
 * @file linalg.cpp
 * @brief Implementation of dense vector linear algebra wrappers over Eigen::Map.
 */
#include "core/algebra/linalg.hpp"

#include <Eigen/Core>
#include <stdexcept>

namespace lsquant {
namespace linalg {

void axpy(Scalar a, const Vector& x, Vector& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("linalg::axpy: vector size mismatch");
  }
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  Eigen::Map<Eigen::Vector<Scalar, -1>> eig_y(y.data(), static_cast<Eigen::Index>(y.size()));
  eig_y = a * eig_x + eig_y;
}

void scal(Scalar a, Vector& x) {
  Eigen::Map<Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  eig_x *= a;
}

void copy(const Vector& x, Vector& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("linalg::copy: vector size mismatch");
  }
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  Eigen::Map<Eigen::Vector<Scalar, -1>> eig_y(y.data(), static_cast<Eigen::Index>(y.size()));
  eig_y = eig_x;
}

Scalar dot(const Vector& x, const Vector& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("linalg::dot: vector size mismatch");
  }
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_y(y.data(), static_cast<Eigen::Index>(y.size()));
  return eig_x.dot(eig_y);
}

double norm2(const Vector& x) {
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  return eig_x.squaredNorm();
}

}  // namespace linalg
}  // namespace lsquant
