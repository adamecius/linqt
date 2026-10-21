/**
 * @file block_diagonal_operator.cpp
 * @brief Implementation of BlockDiagonalOperator.
 */
#include "core/algebra/block_diagonal_operator.hpp"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>

namespace lsquant {

BlockDiagonalOperator::BlockDiagonalOperator(Index numBlocks, Index blockSize,
                                             std::vector<Scalar> blocks)
    : numBlocks_(numBlocks), blockSize_(blockSize), blocks_(std::move(blocks)) {
  if (numBlocks < 0 || blockSize < 0) {
    throw std::invalid_argument(
        "BlockDiagonalOperator: numBlocks and blockSize must be non-negative");
  }
  const std::size_t expected = static_cast<std::size_t>(numBlocks) *
                               static_cast<std::size_t>(blockSize) *
                               static_cast<std::size_t>(blockSize);
  if (blocks_.size() != expected) {
    throw std::invalid_argument(
        "BlockDiagonalOperator: blocks size must equal numBlocks * blockSize^2");
  }

  const std::size_t w = static_cast<std::size_t>(blockSize_);
  const std::size_t w2 = w * w;
  bool herm = true;
  for (Index ik = 0; ik < numBlocks_ && herm; ++ik) {
    const Scalar* blk = blocks_.data() + static_cast<std::size_t>(ik) * w2;
    for (std::size_t i = 0; i < w; ++i) {
      for (std::size_t j = i; j < w; ++j) {
        const Scalar diff = blk[i * w + j] - std::conj(blk[j * w + i]);
        if (std::abs(diff) > 1e-12) {
          herm = false;
          break;
        }
      }
      if (!herm) {
        break;
      }
    }
  }
  isHermitian_ = herm;
}

Index BlockDiagonalOperator::Dimension() const { return numBlocks_ * blockSize_; }

void BlockDiagonalOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index dim = Dimension();
  if (static_cast<Index>(x.size()) != dim || static_cast<Index>(y.size()) != dim) {
    throw std::invalid_argument("BlockDiagonalOperator::Multiply: vector dimension mismatch");
  }
  if (dim == 0) {
    return;
  }

  const Index nk = numBlocks_;
  const std::size_t w = static_cast<std::size_t>(blockSize_);
  const std::size_t w2 = w * w;
  const Scalar* blkBase = blocks_.data();
  const Scalar* xPtr = x.data();
  Scalar* yPtr = y.data();
  const bool aliased = (xPtr == yPtr);

  if (aliased) {
    #pragma omp parallel for schedule(static)
    for (Index ik = 0; ik < nk; ++ik) {
      const std::size_t offset = static_cast<std::size_t>(ik) * w;
      const Scalar* blk = blkBase + static_cast<std::size_t>(ik) * w2;
      Vector xLocal(xPtr + offset, xPtr + offset + w);
      Scalar* yk = yPtr + offset;
      for (std::size_t alpha = 0; alpha < w; ++alpha) {
        const Scalar* row = blk + alpha * w;
        Scalar acc(0.0, 0.0);
        for (std::size_t beta = 0; beta < w; ++beta) {
          acc += row[beta] * xLocal[beta];
        }
        yk[alpha] = a * acc + b * yk[alpha];
      }
    }
  } else {
    #pragma omp parallel for schedule(static)
    for (Index ik = 0; ik < nk; ++ik) {
      const std::size_t offset = static_cast<std::size_t>(ik) * w;
      const Scalar* blk = blkBase + static_cast<std::size_t>(ik) * w2;
      const Scalar* xk = xPtr + offset;
      Scalar* yk = yPtr + offset;
      for (std::size_t alpha = 0; alpha < w; ++alpha) {
        const Scalar* row = blk + alpha * w;
        Scalar acc(0.0, 0.0);
        for (std::size_t beta = 0; beta < w; ++beta) {
          acc += row[beta] * xk[beta];
        }
        yk[alpha] = a * acc + b * yk[alpha];
      }
    }
  }
}

BlockDiagonalOperator BlockDiagonalOperator::Adjoint() const {
  const std::size_t nk = static_cast<std::size_t>(numBlocks_);
  const std::size_t w = static_cast<std::size_t>(blockSize_);
  const std::size_t w2 = w * w;
  std::vector<Scalar> adj(blocks_.size(), Scalar(0.0, 0.0));
  for (std::size_t ik = 0; ik < nk; ++ik) {
    const Scalar* src = blocks_.data() + ik * w2;
    Scalar* dst = adj.data() + ik * w2;
    for (std::size_t i = 0; i < w; ++i) {
      for (std::size_t j = 0; j < w; ++j) {
        dst[i * w + j] = std::conj(src[j * w + i]);
      }
    }
  }
  return BlockDiagonalOperator(numBlocks_, blockSize_, std::move(adj));
}

}  // namespace lsquant
