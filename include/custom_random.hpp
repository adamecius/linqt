#ifndef CUSTOM_RANDOM_HPP
#define CUSTOM_RANDOM_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

namespace custom_random
{
	typedef boost::random::uniform_real_distribution<my::real> uniform_real_dist;
	typedef boost::random::mt19937 generator;
};

static custom_random::generator global_generator;

#endif
