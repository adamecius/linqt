/**
 * @file contraction.hpp
 * @brief Rank-1 contraction of Chebyshev moment vectors against coefficient tables.
 *
 * Evaluates the point-by-point series contraction `values[i] = sum_m table[m, i] * moments[m]`
 * with a deterministic sequential sum over the moment index `m` for each evaluation point `i`,
 * parallelised across `i` with OpenMP static scheduling.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_CONTRACTION_HPP_
#define LSQUANT_CORE_CHEBYSHEV_CONTRACTION_HPP_

#include <vector>

#include "core/algebra/types.hpp"
#include "core/chebyshev/chebyshev_moments.hpp"

namespace lsquant {

/**
 * @brief Contract a rank-1 real moment vector with a 2D coefficient table into a buffer.
 *
 * Solves the problem of evaluating a 1D Chebyshev series at `numPoints` grid coordinates
 * from precomputed expansion coefficients `table[m, i]` and damped real moments `moments[m]`
 * while preserving the exact sequential IEEE-754 summation order over `m = 0 .. numMoments - 1`
 * required for bit-for-bit reproducibility across all OpenMP thread counts. Both `moments` and
 * `table` are dimensionless in the rescaled Chebyshev domain `[-1, 1]`; physical prefactors are
 * applied by the caller. Complexity is `O(numMoments * numPoints)` floating-point operations and
 * `O(1)` auxiliary memory.
 *
 * @param moments      Pointer to the first element of the real moment vector of length `numMoments`.
 * @param momentStride Element stride between consecutive moments in `moments` (e.g., `1` for a
 *                     contiguous `double` array, `2` for the real view of a complex array).
 * @param numMoments   Number of expansion terms `M >= 0`.
 * @param table        Pointer to the first element of the `numMoments x numPoints` coefficient table.
 * @param tableStrideM Element stride along the moment axis `m` of `table` (`numPoints` for row-major).
 * @param tableStrideN Element stride along the grid axis `i` of `table` (`1` for row-major).
 * @param numPoints    Number of evaluation grid points `N >= 0`.
 * @param values       Output buffer of length `numPoints` receiving `sum_m table[m, i] * moments[m]`.
 *
 * @warning Throws `std::invalid_argument` if `numMoments < 0`, `numPoints < 0`, or any required
 *          pointer is null when its corresponding extent is positive.
 * @par Pitfall
 *      Never parallelise or vector-reduce the inner loop over `m`; reordering the `m` summation
 *      changes IEEE-754 rounding in the 6th significant digit of formatted output curves.
 */
void ContractRank1(const double* moments, Index momentStride, Index numMoments,
                   const double* table, Index tableStrideM, Index tableStrideN,
                   Index numPoints, double* values);

/**
 * @brief Contract a contiguous rank-1 real moment vector with a row-major `(M, N)` table.
 *
 * Evaluates `values[i] = sum_{m=0}^{M-1} table[m * N + i] * moments[m]` sequentially in `m`
 * for each grid index `i in 0 .. N - 1`. Complexity is `O(M * N)` time and `O(N)` output space.
 *
 * @param momentsContiguous real moment vector of length `M`.
 * @param table   Row-major coefficient table of length `M * N`.
 * @return        Vector of `N` contracted values.
 *
 * @warning Throws `std::invalid_argument` if `moments` is empty when `table` is non-empty or if
 *          `table.size()` is not an integer multiple of `moments.size()`.
 * @par Pitfall
 *      `table` must be stored in row-major `(M, N)` layout where the moment index `m` is the
 *      slow axis and the evaluation point index `i` is the fast axis.
 */
std::vector<double> ContractRank1(const std::vector<double>& moments,
                                  const std::vector<double>& table);

/**
 * @brief Contract the real part of a rank-1 `ChebyshevMoments` tensor with a row-major `(M, N)` table.
 *
 * Extracts the real components of `moments.Data()` (with stride 2 over the underlying complex
 * storage) and contracts them against the `M x N` row-major `table`. Complexity is `O(M * N)`.
 *
 * @param moments Rank-1 `ChebyshevMoments` tensor of length `M`.
 * @param table   Row-major coefficient table of length `M * N`.
 * @return        Vector of `N` contracted values.
 *
 * @warning Throws `std::invalid_argument` if `moments.Rank() != 1` or if `table.size()` is not
 *          divisible by `moments.Dims()[0]`.
 * @par Pitfall
 *      Uses `moments` as stored; any kernel damping or truncation must be applied before calling
 *      this overload.
 */
std::vector<double> ContractRank1(const ChebyshevMoments& moments,
                                  const std::vector<double>& table);

namespace chebyshev {
using lsquant::ContractRank1;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_CONTRACTION_HPP_
