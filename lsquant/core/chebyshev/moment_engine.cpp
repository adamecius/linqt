/**
 * @file moment_engine.cpp
 * @brief Implementation of MomentEngine::Run.
 */
#include "core/chebyshev/moment_engine.hpp"

#include <algorithm>
#include <stdexcept>
#include <variant>
#include <vector>

#include "core/algebra/linalg.hpp"
#include "core/chebyshev/chebyshev_recursion.hpp"
#include "core/chebyshev/chebyshev_vectors.hpp"

namespace lsquant {

namespace {

struct SideProgram {
  std::vector<int> preOperators;
  bool hasExpand = false;
  int numMoments = 0;
  std::vector<int> postOperators;
};

void ApplySingleVectorOperators(const std::vector<const LinearOperator*>& operators,
                                const std::vector<int>& opIndices, const Vector& init,
                                Vector& out, Vector& scratch) {
  linalg::copy(init, out);
  for (int idx : opIndices) {
    const LinearOperator& op = *operators[static_cast<std::size_t>(idx)];
    if (op.IsIdentity()) {
      continue;
    }
    linalg::copy(out, scratch);
    op.Multiply(Scalar(1.0, 0.0), scratch, Scalar(0.0, 0.0), out);
  }
}

void AccumulateWeightedDot(Scalar& target, const Scalar& weight, const Vector& bra,
                           const Vector& ket) {
  const Scalar d = linalg::dot(bra, ket);
  if (weight.imag() == 0.0) {
    target += weight.real() * d;
  } else {
    target += weight * d;
  }
}

int PaddedMoments(int numMoments, int batchSize) {
  return ((numMoments + batchSize - 1) / batchSize) * batchSize;
}

}  // namespace

MomentEngine::MomentEngine(const LinearOperator& rescaledH,
                           std::vector<const LinearOperator*> operators,
                           StateGenerator& states, int batchSize)
    : rescaledH_(&rescaledH),
      operators_(std::move(operators)),
      states_(&states),
      batchSize_(batchSize) {
  if (rescaledH.Dimension() <= 0) {
    throw std::invalid_argument("MomentEngine: rescaledH dimension must be positive");
  }
  if (states.Dimension() != rescaledH.Dimension()) {
    throw std::invalid_argument("MomentEngine: states dimension mismatch");
  }
  if (batchSize <= 0) {
    throw std::invalid_argument("MomentEngine: batchSize must be positive");
  }
  for (const LinearOperator* op : operators_) {
    if (op == nullptr) {
      throw std::invalid_argument("MomentEngine: null operator pointer");
    }
    if (op->Dimension() != rescaledH.Dimension()) {
      throw std::invalid_argument("MomentEngine: operator dimension mismatch");
    }
  }
}

ChebyshevMoments MomentEngine::Run(const Plan& plan) {
  if (plan.code.empty()) {
    throw std::invalid_argument("MomentEngine::Run: plan code must not be empty");
  }

  SideProgram ketProg;
  SideProgram braProg;
  SideProgram* activeSide = &ketProg;
  bool seenMirror = false;
  bool seenTrace = false;
  std::vector<Scalar> traceWeights;

  for (std::size_t pc = 0; pc < plan.code.size(); ++pc) {
    if (seenTrace) {
      throw std::invalid_argument("MomentEngine::Run: Trace must be the final instruction");
    }
    const Instruction& inst = plan.code[pc];
    if (std::holds_alternative<ApplySeries>(inst)) {
      throw std::runtime_error("ApplySeries: not implemented in stage 1");
    } else if (const auto* app = std::get_if<ApplyOperator>(&inst)) {
      if (app->operatorIndex < 0 ||
          app->operatorIndex >= static_cast<int>(operators_.size())) {
        throw std::invalid_argument("MomentEngine::Run: operatorIndex out of range");
      }
      if (!activeSide->hasExpand) {
        activeSide->preOperators.push_back(app->operatorIndex);
      } else {
        activeSide->postOperators.push_back(app->operatorIndex);
      }
    } else if (const auto* exp = std::get_if<Expand>(&inst)) {
      if (exp->numMoments <= 0) {
        throw std::invalid_argument("MomentEngine::Run: Expand numMoments must be positive");
      }
      if (activeSide->hasExpand) {
        throw std::runtime_error("Expand: rank > 2 not implemented in stage 1");
      }
      activeSide->hasExpand = true;
      activeSide->numMoments = exp->numMoments;
    } else if (std::holds_alternative<MirrorToBra>(inst)) {
      if (seenMirror) {
        throw std::invalid_argument("MomentEngine::Run: MirrorToBra may appear at most once");
      }
      seenMirror = true;
      activeSide = &braProg;
    } else if (const auto* tr = std::get_if<Trace>(&inst)) {
      seenTrace = true;
      traceWeights = tr->weights;
    }
  }

  if (!seenTrace) {
    throw std::invalid_argument("MomentEngine::Run: plan must end with a Trace instruction");
  }

  std::vector<int> dims;
  if (braProg.hasExpand) {
    dims.push_back(braProg.numMoments);
  }
  if (ketProg.hasExpand) {
    dims.push_back(ketProg.numMoments);
  }

  ChebyshevMoments mu(dims);
  mu.systemSize = rescaledH_->Dimension();
  mu.numberOfStates = static_cast<int>(states_->NumberOfStates());
  mu.stateLabel = states_->Label();
  mu.word = plan.word;

  const std::size_t tensorSize = mu.Data().size();
  if (traceWeights.size() != 1 && traceWeights.size() != tensorSize) {
    throw std::invalid_argument(
        "MomentEngine::Run: Trace weights size must be 1 or match tensor size");
  }
  const bool scalarWeight = (traceWeights.size() == 1);

  const Index dim = rescaledH_->Dimension();
  states_->Reset();

  Vector r(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector ketInit(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector braInit(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector scratch(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));

  if (!braProg.hasExpand && !ketProg.hasExpand) {
    // Rank 0: neither side opens a Chebyshev index.
    while (states_->Next(r)) {
      ApplySingleVectorOperators(operators_, ketProg.preOperators, r, ketInit, scratch);
      ApplySingleVectorOperators(operators_, braProg.preOperators, r, braInit, scratch);
      AccumulateWeightedDot(mu.Data()[0], traceWeights[0], braInit, ketInit);
    }
    return mu;
  }

  if (!braProg.hasExpand && ketProg.hasExpand) {
    // Rank 1 (ket expanded): matches legacy SpectralMoments when postOperators is empty.
    const int numMoms = ketProg.numMoments;
    if (ketProg.postOperators.empty()) {
      ChebyshevRecursion rec(*rescaledH_);
      while (states_->Next(r)) {
        ApplySingleVectorOperators(operators_, braProg.preOperators, r, braInit, scratch);
        ApplySingleVectorOperators(operators_, ketProg.preOperators, r, ketInit, scratch);
        rec.Start(ketInit);
        for (int m = 0; m < numMoms; ++m) {
          const Scalar& w = scalarWeight ? traceWeights[0]
                                         : traceWeights[static_cast<std::size_t>(m)];
          AccumulateWeightedDot(mu.At(m), w, braInit, rec.Current());
          rec.Step();
        }
      }
    } else {
      const int effBatch = std::min(batchSize_, numMoms);
      const int paddedM = PaddedMoments(numMoms, effBatch);
      ChebyshevRecursion rec(*rescaledH_);
      ChebyshevVectors vecR(dim, effBatch);
      while (states_->Next(r)) {
        ApplySingleVectorOperators(operators_, braProg.preOperators, r, braInit, scratch);
        ApplySingleVectorOperators(operators_, ketProg.preOperators, r, ketInit, scratch);
        rec.Start(ketInit);
        for (int mR = 0; mR < paddedM; mR += effBatch) {
          vecR.Fill(rec);
          for (int opIdx : ketProg.postOperators) {
            vecR.ApplyToAll(*operators_[static_cast<std::size_t>(opIdx)]);
          }
          for (int m0 = 0; m0 < effBatch; ++m0) {
            const int m = mR + m0;
            if (m < numMoms) {
              const Scalar& w = scalarWeight ? traceWeights[0]
                                             : traceWeights[static_cast<std::size_t>(m)];
              AccumulateWeightedDot(mu.At(m), w, braInit, vecR.At(m0));
            }
          }
        }
      }
    }
    return mu;
  }

  if (braProg.hasExpand && !ketProg.hasExpand) {
    // Rank 1 (bra expanded).
    const int numMoms = braProg.numMoments;
    const int effBatch = std::min(batchSize_, numMoms);
    const int paddedM = PaddedMoments(numMoms, effBatch);
    ChebyshevRecursion recL(*rescaledH_);
    ChebyshevVectors vecL(dim, effBatch);
    while (states_->Next(r)) {
      ApplySingleVectorOperators(operators_, ketProg.preOperators, r, ketInit, scratch);
      ApplySingleVectorOperators(operators_, braProg.preOperators, r, braInit, scratch);
      recL.Start(braInit);
      for (int mL = 0; mL < paddedM; mL += effBatch) {
        vecL.Fill(recL);
        for (int opIdx : braProg.postOperators) {
          vecL.ApplyToAll(*operators_[static_cast<std::size_t>(opIdx)]);
        }
        for (int m0 = 0; m0 < effBatch; ++m0) {
          const int m = mL + m0;
          if (m < numMoms) {
            const Scalar& w = scalarWeight ? traceWeights[0]
                                           : traceWeights[static_cast<std::size_t>(m)];
            AccumulateWeightedDot(mu.At(m), w, vecL.At(m0), ketInit);
          }
        }
      }
    }
    return mu;
  }

  // Rank 2 (both bra and ket expanded): matches legacy CorrelationExpansionMoments block order.
  const int numMomsL = braProg.numMoments;
  const int numMomsR = ketProg.numMoments;
  const int effBatchL = std::min(batchSize_, numMomsL);
  const int effBatchR = std::min(batchSize_, numMomsR);
  const int paddedML = PaddedMoments(numMomsL, effBatchL);
  const int paddedMR = PaddedMoments(numMomsR, effBatchR);

  ChebyshevRecursion recL(*rescaledH_);
  ChebyshevRecursion recR(*rescaledH_);
  ChebyshevVectors vecL(dim, effBatchL);
  ChebyshevVectors vecR(dim, effBatchR);

  while (states_->Next(r)) {
    ApplySingleVectorOperators(operators_, braProg.preOperators, r, braInit, scratch);
    recL.Start(braInit);
    for (int mL = 0; mL < paddedML; mL += effBatchL) {
      vecL.Fill(recL);
      for (int opIdx : braProg.postOperators) {
        vecL.ApplyToAll(*operators_[static_cast<std::size_t>(opIdx)]);
      }

      ApplySingleVectorOperators(operators_, ketProg.preOperators, r, ketInit, scratch);
      recR.Start(ketInit);
      for (int mR = 0; mR < paddedMR; mR += effBatchR) {
        vecR.Fill(recR);
        for (int opIdx : ketProg.postOperators) {
          vecR.ApplyToAll(*operators_[static_cast<std::size_t>(opIdx)]);
        }

        for (int m0 = 0; m0 < effBatchL; ++m0) {
          const int row = mL + m0;
          if (row >= numMomsL) {
            continue;
          }
          for (int m1 = 0; m1 < effBatchR; ++m1) {
            const int col = mR + m1;
            if (col >= numMomsR) {
              continue;
            }
            const std::size_t wIdx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(numMomsR) +
                static_cast<std::size_t>(col);
            const Scalar& w = scalarWeight ? traceWeights[0] : traceWeights[wIdx];
            AccumulateWeightedDot(mu.At(row, col), w, vecL.At(m0), vecR.At(m1));
          }
        }
      }
    }
  }

  return mu;
}

}  // namespace lsquant
