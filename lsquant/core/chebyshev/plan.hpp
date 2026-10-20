/**
 * @file plan.hpp
 * @brief Five-instruction bytecode representation of a KPM operator word.
 *
 * Defines the physical-name-free instruction set (`ApplyOperator`, `Expand`,
 * `ApplySeries`, `MirrorToBra`, `Trace`) and `Plan` container executed by
 * `MomentEngine::Run`.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_PLAN_HPP_
#define LSQUANT_CORE_CHEBYSHEV_PLAN_HPP_

#include <string>
#include <variant>
#include <vector>

#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Instruction applying registered static operator `O_{operatorIndex}` to every live vector.
 *
 * Replaces each live vector `v` on the active side (ket or bra) with `O_k v`.
 *
 * @warning `operatorIndex` must be in `[0, operators.size())` of the executing `MomentEngine`.
 * @par Pitfall
 *      When applied after `Expand`, transforms every vector in the active Chebyshev batch.
 */
struct ApplyOperator {
  int operatorIndex;
};

/**
 * @brief Instruction opening one tensor index of length `numMoments` via Chebyshev recursion.
 *
 * Replaces a live vector `v` on the active side with the sequence `{T_0(H~) v, ..., T_{M-1}(H~) v}`.
 *
 * @warning `numMoments` must be strictly positive; at most one `Expand` per side is supported
 *          in stage 1.
 * @par Pitfall
 *      Each `Expand` increases the rank of the output `ChebyshevMoments` tensor by 1.
 */
struct Expand {
  int numMoments;
};

/**
 * @brief Instruction contracting a Chebyshev series `sum_m coeffs[m] T_m(H~) v` on the fly.
 *
 * Present in the `Instruction` variant type to fix the five-instruction bytecode grammar;
 * execution belongs to stage 4 (`MomentEngine::Run` throws in stage 1).
 *
 * @warning Executing `ApplySeries` in stage 1 throws `std::runtime_error("ApplySeries: not implemented in stage 1")`.
 * @par Pitfall
 *      Do not emit `ApplySeries` in stage-1 plans.
 */
struct ApplySeries {
  std::vector<Scalar> coeffs;
};

/**
 * @brief Instruction switching execution from the ket side to the bra side on a fresh copy of `|r>`.
 *
 * May appear at most once in a `Plan`. All instructions prior to `MirrorToBra` act on the ket
 * side; instructions between `MirrorToBra` and `Trace` act on the bra side starting from `|r>`.
 *
 * @warning Throws `std::invalid_argument` if `MirrorToBra` appears more than once in a `Plan`.
 * @par Pitfall
 *      If `MirrorToBra` is omitted, the bra side remains the unmodified initial trace state `|r>`.
 */
struct MirrorToBra {};

/**
 * @brief Terminal instruction contracting live bra and ket vectors into the moment tensor.
 *
 * Accumulates `mu[i_bra..., i_ket...] += weights[...] * dot(bra_{i_bra}, ket_{i_ket})` for each
 * trace state. `weights` must have size `1` (uniform scalar weight) or match the total element
 * count of `mu`.
 *
 * @warning Must appear as the final instruction of every `Plan`.
 * @par Pitfall
 *      Real weights are multiplied as `double * Scalar` to preserve IEEE-754 signed zeros in
 *      imaginary components.
 */
struct Trace {
  std::vector<Scalar> weights;
};

/**
 * @brief Variant holding any one of the five engine instructions.
 */
using Instruction = std::variant<ApplyOperator, Expand, ApplySeries, MirrorToBra, Trace>;

/**
 * @brief Compiled instruction sequence and descriptive word string for `MomentEngine::Run`.
 *
 * Pairs the executable `code` vector with the human-readable `word` label stored in the
 * resulting `ChebyshevMoments` metadata.
 *
 * @warning `code` must end with a single `Trace` instruction.
 * @par Pitfall
 *      Operator indices in `code` refer positionally to the `operators` vector passed to
 *      `MomentEngine`.
 */
struct Plan {
  std::vector<Instruction> code;
  std::string word;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_PLAN_HPP_
