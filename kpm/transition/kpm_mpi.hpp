/*
 * kpm_mpi.hpp
 *
 *  Created on: 24/08/2016
 *      Author: jgarcia
 */

#ifndef KPM_MPI_HPP_
#define KPM_MPI_HPP_

/*
 * kpm_mpi_macros.hpp
 *
 *  Created on: May 27, 2016
 *      Author: jgarcia
 */


#include "mpi.h"


#define KPM_MPI_INIT() MPI::Init ()
#define KPM_MPI_FINALIZE()   MPI::Finalize ( )
#define KPM_MPI_GETRANK() MPI::COMM_WORLD.Get_rank ( )
#define KPM_MPI_GETPROC()  MPI::COMM_WORLD.Get_size ( )


#endif /* KPM_KPM_MPI_HPP_ */
