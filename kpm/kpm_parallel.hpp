/*
 * kpm_parallel.hpp
 *
 *  Created on: 24/08/2016
 *      Author: jgarcia
 */

#ifndef KPM_KPM_PARALLEL_HPP_
#define KPM_KPM_PARALLEL_HPP_


#include "kpm_mpi.hpp"

#include "kpm_openmp.hpp"

#include "kpm_cuda.hpp"

#define KPM_PARALLEL_INIT() MPI::Init(); omp_init();

#define KPM_PARALLEL_FINALIZE()  MPI::Finalize ();

#endif /* KPM_KPM_PARALLEL_HPP_ */
