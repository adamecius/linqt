/**
 * @file state_generator.hpp
 * @brief Source-state generators for exact and stochastic trace evaluations.
 *
 * Defines the polymorphic `StateGenerator` interface and its concrete implementations:
 * `RandomPhaseStates` (reproducing the legacy Mersenne Twister random-phase convention),
 * `LocalStates` (unit coordinate vectors for local projections or full exact traces),
 * and `UserStates` (caller-provided state vectors).
 */
#ifndef LSQUANT_CORE_ALGEBRA_STATE_GENERATOR_HPP_
#define LSQUANT_CORE_ALGEBRA_STATE_GENERATOR_HPP_

#include <random>
#include <string>
#include <vector>

#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Abstract iterator-style generator yielding quantum state vectors for trace estimation.
 *
 * Solves the problem of streaming source states `|r>` one at a time into the moment engine
 * without storing all `NumberOfStates()` vectors in memory simultaneously. Each call to
 * `Next(out)` populates `out` with a vector of length `Dimension()` and advances the cursor;
 * `Reset()` rewinds the generator to the initial state so multiple plans can be evaluated over
 * the identical state sequence.
 *
 * @warning `Next` mutates the generator's internal state and is not meant to be called
 *          concurrently from multiple threads.
 * @par Pitfall
 *      Always call `Reset()` before reusing a `StateGenerator` across successive `MomentEngine`
 *      runs so that the cursor starts at index 0.
 */
class StateGenerator {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~StateGenerator() = default;

  /**
   * @brief Return the Hilbert-space dimension `N` of generated vectors.
   *
   * @return Vector length `N > 0`.
   */
  virtual Index Dimension() const = 0;

  /**
   * @brief Return the total number of state vectors in the sequence.
   *
   * @return Count of states `R >= 0`.
   */
  virtual Index NumberOfStates() const = 0;

  /**
   * @brief Populate `out` with the next state vector in the sequence.
   *
   * @param out Output vector resized/filled with `Dimension()` complex entries.
   * @return    `true` if a state was produced, `false` when the sequence is exhausted.
   */
  virtual bool Next(Vector& out) = 0;

  /**
   * @brief Rewind the generator to its initial state and re-seed any internal PRNG.
   */
  virtual void Reset() = 0;

  /**
   * @brief Return the canonical state-set label used in output filenames and headers.
   *
   * @return One of `"default"`, `"exact"`, `"local"`, or `"user"`.
   */
  virtual std::string Label() const = 0;
};

/**
 * @brief Stochastic-trace generator yielding normalised random-phase vectors.
 *
 * Reproduces the exact legacy random-phase convention (`RANDOM_STATE`):
 * for each component `i in 0..N-1`, draws `u` from `std::uniform_real_distribution<double>(0.0, 1.0)`
 * driven by `std::mt19937(seed)`, sets `phi = 2 * pi * u`, and assigns
 * `x_i = (cos(phi) / sqrt(N), sin(phi) / sqrt(N))`. Both the engine and the distribution
 * object are constructed once per generator instance. Cost per state is `O(N)`.
 *
 * @warning Changing the distribution construction lifetime or PRNG type changes the random
 *          stream and breaks stochastic-trace regression goldens.
 * @par Pitfall
 *      `Reset()` must both re-seed `rng_` with `seed_` and call `uni_.reset()` so that cached
 *      distribution state cannot leak across runs.
 */
class RandomPhaseStates : public StateGenerator {
public:
  /**
   * @brief Construct a random-phase generator for `numberOfStates` vectors of length `dimension`.
   *
   * @param dimension      Hilbert-space dimension `N > 0`.
   * @param numberOfStates Number of random vectors `R >= 0`.
   * @param seed           Deterministic `std::mt19937` seed (default `12345`).
   *
   * @warning Throws `std::invalid_argument` if `dimension <= 0` or `numberOfStates < 0`.
   * @par Pitfall
   *      Pass the same `seed` to reproduce a previous stochastic trace bit for bit.
   */
  RandomPhaseStates(Index dimension, Index numberOfStates, unsigned seed = 12345);

  /**
   * @brief Return the vector dimension `N`.
   *
   * @return Dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Return the number of random vectors `R`.
   *
   * @return State count `R >= 0`.
   */
  Index NumberOfStates() const override;

  /**
   * @brief Generate the next normalised random-phase vector into `out`.
   *
   * @param out Vector overwritten with the next random-phase state.
   * @return    `true` if generated, `false` if all `NumberOfStates()` have been consumed.
   */
  bool Next(Vector& out) override;

  /**
   * @brief Rewind the counter and reset the `std::mt19937` engine to its initial `seed`.
   */
  void Reset() override;

  /**
   * @brief Return `"default"` (matching legacy `gen.StateLabel()`).
   *
   * @return `"default"`.
   */
  std::string Label() const override;

private:
  Index dimension_;
  Index numberOfStates_;
  unsigned seed_;
  Index current_ = 0;
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uni_;
};

/**
 * @brief Deterministic generator yielding unit coordinate vectors `|e_{s}>`.
 *
 * For each index `s` in `sites`, yields a vector of length `dimension` with `1.0` at
 * index `s` and `0.0` elsewhere. Reports label `"exact"` if and only if `sites` is the
 * complete ascending sequence `0, 1, ..., dimension - 1`, and `"local"` otherwise.
 * Cost per state is `O(dimension)`.
 *
 * @warning Throws `std::invalid_argument` if `dimension <= 0` or any site index is outside
 *          `[0, dimension)`.
 * @par Pitfall
 *      Full exact trace over `N` states requires `O(N)` polynomial recursions; use only on
 *      small or medium systems.
 */
class LocalStates : public StateGenerator {
public:
  /**
   * @brief Construct a local-state generator for the specified site indices.
   *
   * @param dimension Hilbert-space dimension `N > 0`.
   * @param sites     Ordered list of site indices in `[0, dimension)`.
   *
   * @warning Throws `std::invalid_argument` if `dimension <= 0` or any site is out of bounds.
   * @par Pitfall
   *      The label is `"exact"` only when `sites` equals `0..dimension-1` in ascending order.
   */
  LocalStates(Index dimension, std::vector<Index> sites);

  /**
   * @brief Return the vector dimension `N`.
   *
   * @return Dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Return `sites.size()`.
   *
   * @return Number of local site states.
   */
  Index NumberOfStates() const override;

  /**
   * @brief Populate `out` with the next unit vector `|e_{sites[current]}>`.
   *
   * @param out Vector overwritten with zeros and `1.0` at `sites[current]`.
   * @return    `true` if a state was written, `false` when `sites` is exhausted.
   */
  bool Next(Vector& out) override;

  /**
   * @brief Rewind the site cursor to `0`.
   */
  void Reset() override;

  /**
   * @brief Return `"exact"` iff `sites == 0..N-1`, else `"local"`.
   *
   * @return `"exact"` or `"local"`.
   */
  std::string Label() const override;

private:
  Index dimension_;
  std::vector<Index> sites_;
  Index current_ = 0;
  std::string label_;
};

/**
 * @brief Generator yielding caller-supplied state vectors in order.
 *
 * Stores an explicit list of dense `Vector` states of uniform positive dimension and returns
 * them sequentially on `Next(out)`. Reports label `"user"`. Cost per state is `O(dimension)`.
 *
 * @warning Throws `std::invalid_argument` if the vector list is empty, if vectors have length
 *          `0`, or if vector lengths are inconsistent.
 * @par Pitfall
 *      Caller vectors are not automatically normalised; supply normalised states if unit norm
 *      is expected.
 */
class UserStates : public StateGenerator {
public:
  /**
   * @brief Construct a generator over an explicit list of state vectors.
   *
   * @param states Non-empty list of equal-length non-empty `Vector` instances.
   *
   * @warning Throws `std::invalid_argument` if `states` is empty or dimensions mismatch.
   * @par Pitfall
   *      Moves `states` into internal storage; the caller's vector is consumed when passed as
   *      an rvalue.
   */
  explicit UserStates(std::vector<Vector> states);

  /**
   * @brief Return the common length of the stored user vectors.
   *
   * @return Dimension `N > 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Return the number of stored user vectors.
   *
   * @return State count `R > 0`.
   */
  Index NumberOfStates() const override;

  /**
   * @brief Copy the next user vector into `out`.
   *
   * @param out Output vector receiving the next user state.
   * @return    `true` if a state was copied, `false` when all states have been yielded.
   */
  bool Next(Vector& out) override;

  /**
   * @brief Rewind the cursor to the first user vector.
   */
  void Reset() override;

  /**
   * @brief Return `"user"`.
   *
   * @return `"user"`.
   */
  std::string Label() const override;

private:
  Index dimension_ = 0;
  std::vector<Vector> states_;
  Index current_ = 0;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_STATE_GENERATOR_HPP_
