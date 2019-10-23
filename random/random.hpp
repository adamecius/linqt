#ifndef NUMCAL_RANDOM_HPP
#define NUMCAL_RANDOM_HPP

#include "types_definitions.hpp"
#include "types_definitions.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <ctime>
#include <sys/time.h>
#include <unistd.h>

namespace NumCal
{

namespace random
{
	typedef boost::random::uniform_real_distribution< real > uniform_real_dist;
	typedef boost::random::mt19937 generator;


	unsigned int random_seed();


};

};

#endif
