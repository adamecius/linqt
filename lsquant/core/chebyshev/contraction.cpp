/**
 * @file contraction.cpp
 * @brief Implementation of deterministic rank-1 Chebyshev series contraction.
 */
#include "core/chebyshev/contraction.hpp"

#include <stdexcept>

namespace lsquant {

void ContractRank1(const double* moments, Index momentStride, Index numMoments,
                   const double* table, Index tableStrideM, Index tableStrideN,
                   Index numPoints, double* values) {
  if (numMoments < 0) {
    throw std::invalid_argument("ContractRank1: numMoments must be non-negative");
  }
  if (numPoints < 0) {
    throw std::invalid_argument("ContractRank1: numPoints must be non-negative");
  }
  if (numMoments > 0 && moments == nullptr) {
    throw std::invalid_argument("ContractRank1: moments pointer is null");
  }
  if (numMoments > 0 && numPoints > 0 && table == nullptr) {
    throw std::invalid_argument("ContractRank1: table pointer is null");
  }
  if (numPoints > 0 && values == nullptr) {
    throw std::invalid_argument("ContractRank1: values pointer is null");
  }

#pragma omp parallel for schedule(static)
  for (Index i = 0; i < numPoints; ++i) {
    double s = 0.0;
    const double* colPtr = table + i * tableStrideN;
    for (Index m = 0; m < numMoments; ++m) {
      s += colPtr[m * tableStrideM] * moments[m * momentStride];
    }
    values[i] = s;
  }
}

std::vector<double> ContractRank1(const std::vector<double>& moments,
                                  const std::vector<double>& table) {
  const Index numMoments = static_cast<Index>(moments.size());
  if (numMoments == 0) {
    if (!table.empty()) {
      throw std::invalid_argument(
          "ContractRank1: table must be empty when moments vector is empty");
    }
    return {};
  }
  const Index totalSize = static_cast<Index>(table.size());
  if (totalSize % numMoments != 0) {
    throw std::invalid_argument(
        "ContractRank1: table size must be a multiple of moments size");
  }
  const Index numPoints = totalSize / numMoments;
  std::vector<double> values(static_cast<std::size_t>(numPoints), 0.0);
  ContractRank1(moments.data(), 1, numMoments, table.data(), numPoints, 1, numPoints,
                values.data());
  return values;
}

std::vector<double> ContractRank1(const ChebyshevMoments& moments,
                                  const std::vector<double>& table) {
  if (moments.Rank() != 1) {
    throw std::invalid_argument("ContractRank1: ChebyshevMoments must have rank 1");
  }
  const Index numMoments = static_cast<Index>(moments.Dims()[0]);
  const Index totalSize = static_cast<Index>(table.size());
  if (numMoments <= 0 || totalSize % numMoments != 0) {
    throw std::invalid_argument(
        "ContractRank1: table size must be a multiple of moments dimension");
  }
  const Index numPoints = totalSize / numMoments;
  std::vector<double> values(static_cast<std::size_t>(numPoints), 0.0);
  const auto* realPtr = reinterpret_cast<const double*>(moments.Data().data());
  ContractRank1(realPtr, 2, numMoments, table.data(), numPoints, 1, numPoints, values.data());
  return values;
}

}  // namespace lsquant
