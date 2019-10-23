/*
 * kpm_openmp.hpp
 *
 *  Created on: 24/08/2016
 *      Author: jgarcia
 */

#ifndef KPM_OPENMP_HPP_
#define KPM_OPENMP_HPP_



//the compiler defines the macro _OPENMP when the openmp option is used
#include <omp.h>
#include <cstdlib>

const unsigned int
KPM_OMP_NUM_THREADS = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
#define KPM_TIMMING omp_get_wtime

#endif /* KPM_KPM_OPENMP_HPP_ */
