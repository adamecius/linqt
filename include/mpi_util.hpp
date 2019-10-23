/*
 * kpm_mpi.hpp
 *
 *  Created on: 24/08/2016
 *      Author: jgarcia
 */

#ifndef MPI_UTIL_HPP_
#define MPI_UTIL_HPP_

/*
 * kpm_mpi_macros.hpp
 *
 *  Created on: May 27, 2016
 *      Author: jgarcia
 */

#include "mpi.h"
#include <cmath>
#include <iostream>
#include <string>

#define NUMCAL_MPI_INIT() MPI::Init ()
#define NUMCAL_MPI_FINALIZE()   MPI::Finalize ( )
#define NUMCAL_MPI_GETRANK() MPI::COMM_WORLD.Get_rank ( )
#define NUMCAL_MPI_GETPROC()  MPI::COMM_WORLD.Get_size ( )

//#define NUM_CAL_SEND_INT(a) 	MPI::Comm::Bcast(&a,1, MPI::INT , 0, MPI::COMM_WORLD)
namespace NumCal
{

class cout_class
{

public:
	cout_class(){ }

	cout_class(const std::string _init_message)
	{
		if( NUMCAL_MPI_GETRANK() == 0)
			std::cout << _init_message ;
	};

	template <typename T>
	cout_class &operator << (const T &v)
	{
		if( NUMCAL_MPI_GETRANK() == 0)
			std::cout << v;
		return *this;
	};

};

class cerr_class
{

public:
	cerr_class(){ }

	cerr_class(const std::string _init_message)
	{
		if( NUMCAL_MPI_GETRANK() == 0)
			std::cerr << _init_message ;
	};

	template <typename T>
	cerr_class &operator << (const T &v)
	{
		if( NUMCAL_MPI_GETRANK() == 0)
			std::cerr << v;
		return *this;
	};

};
static cout_class cout;
static cerr_class cerr;
static std::string endl("\n");


}
#endif /* KPM_KPM_MPI_HPP_ */
