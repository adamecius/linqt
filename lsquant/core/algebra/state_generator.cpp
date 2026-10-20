/**
 * @file state_generator.cpp
 * @brief Implementation of RandomPhaseStates, LocalStates, and UserStates.
 */
#include "core/algebra/state_generator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace lsquant {

RandomPhaseStates::RandomPhaseStates(Index dimension, Index numberOfStates, unsigned seed)
    : dimension_(dimension),
      numberOfStates_(numberOfStates),
      seed_(seed),
      current_(0),
      rng_(static_cast<std::mt19937::result_type>(seed)),
      uni_(0.0, 1.0) {
  if (dimension <= 0) {
    throw std::invalid_argument("RandomPhaseStates: dimension must be positive");
  }
  if (numberOfStates < 0) {
    throw std::invalid_argument("RandomPhaseStates: numberOfStates must be non-negative");
  }
}

Index RandomPhaseStates::Dimension() const { return dimension_; }

Index RandomPhaseStates::NumberOfStates() const { return numberOfStates_; }

bool RandomPhaseStates::Next(Vector& out) {
  if (current_ >= numberOfStates_) {
    return false;
  }
  out.resize(static_cast<std::size_t>(dimension_));
  const double norm = std::sqrt(static_cast<double>(out.size()));
  for (auto& x : out) {
    const double phi = 2.0 * M_PI * uni_(rng_);
    x = Scalar(std::cos(phi) / norm, std::sin(phi) / norm);
  }
  ++current_;
  return true;
}

void RandomPhaseStates::Reset() {
  current_ = 0;
  rng_.seed(static_cast<std::mt19937::result_type>(seed_));
  uni_.reset();
}

std::string RandomPhaseStates::Label() const { return "default"; }

LocalStates::LocalStates(Index dimension, std::vector<Index> sites)
    : dimension_(dimension), sites_(std::move(sites)), current_(0) {
  if (dimension <= 0) {
    throw std::invalid_argument("LocalStates: dimension must be positive");
  }
  for (Index s : sites_) {
    if (s < 0 || s >= dimension) {
      throw std::invalid_argument("LocalStates: site index out of range");
    }
  }
  bool isExact = (static_cast<Index>(sites_.size()) == dimension);
  if (isExact) {
    for (Index i = 0; i < dimension; ++i) {
      if (sites_[static_cast<std::size_t>(i)] != i) {
        isExact = false;
        break;
      }
    }
  }
  label_ = isExact ? "exact" : "local";
}

Index LocalStates::Dimension() const { return dimension_; }

Index LocalStates::NumberOfStates() const { return static_cast<Index>(sites_.size()); }

bool LocalStates::Next(Vector& out) {
  if (current_ >= static_cast<Index>(sites_.size())) {
    return false;
  }
  out.assign(static_cast<std::size_t>(dimension_), Scalar(0.0, 0.0));
  const Index site = sites_[static_cast<std::size_t>(current_)];
  out[static_cast<std::size_t>(site)] = Scalar(1.0, 0.0);
  ++current_;
  return true;
}

void LocalStates::Reset() { current_ = 0; }

std::string LocalStates::Label() const { return label_; }

UserStates::UserStates(std::vector<Vector> states)
    : states_(std::move(states)), current_(0) {
  if (states_.empty()) {
    throw std::invalid_argument("UserStates: state list must not be empty");
  }
  dimension_ = static_cast<Index>(states_[0].size());
  if (dimension_ <= 0) {
    throw std::invalid_argument("UserStates: state vector dimension must be positive");
  }
  for (const auto& v : states_) {
    if (static_cast<Index>(v.size()) != dimension_) {
      throw std::invalid_argument("UserStates: inconsistent state vector dimensions");
    }
  }
}

Index UserStates::Dimension() const { return dimension_; }

Index UserStates::NumberOfStates() const { return static_cast<Index>(states_.size()); }

bool UserStates::Next(Vector& out) {
  if (current_ >= static_cast<Index>(states_.size())) {
    return false;
  }
  out = states_[static_cast<std::size_t>(current_)];
  ++current_;
  return true;
}

void UserStates::Reset() { current_ = 0; }

std::string UserStates::Label() const { return "user"; }

}  // namespace lsquant
