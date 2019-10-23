/*
 * kpm_utilities.hpp
 *
 *  Created on: 25/08/2016
 *      Author: jgarcia
 */

#ifndef KPM_KPM_UTILITIES_HPP_
#define KPM_KPM_UTILITIES_HPP_

// Defines the types of the program
#include "types_definitions.hpp"
// Defines macros of parallelism
#include "kpm_parallel.hpp"

#include <iostream>
#include <string>


namespace kpm_util
{

my::integer
GetMomentsPerNode( const my::integer numMoms);

void GetMPIStatus(const my::integer numMoms);

void GetOMPStatus();

std::string
GetNodeLabel();

};

#endif /* KPM_KPM_UTILITIES_HPP_ */
