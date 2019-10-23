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


inline
int omp_init()
{
	int num_threads=1;

	if(const char* env_p = std::getenv("OMP_NUM_THREADS"))
		num_threads=std::max(atoi(env_p), 1);
	return  num_threads;
};

const unsigned int
KPM_OMP_NUM_THREADS = omp_init();




#define KPM_TIMMING omp_get_wtime

#endif /* KPM_KPM_OPENMP_HPP_ */
