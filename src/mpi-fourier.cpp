#include <cmath>
#include <cstdlib>     /* srand, rand */
#include <iostream>
#include "lattice_fftw3_mpi.hpp"
#include <gsl/gsl_rng.h>
#include <sys/types.h>
#include <unistd.h>

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// Initialize the lattice fast-fourier mpi-class
	LatticeFFTW PotFFT;
	std::string fftw_label="myfftw";
	PotFFT.Init(MPI_COMM_WORLD, fftw_label);
	
	
	//Get local dimension of the subspace per process 
	int locDim = PotFFT.LocalDim();
	//Get local dimension of the parallelized direction dim0
	int locDim0 = PotFFT.LocalDim0();
	//Get local index initial index of the parallelized direction i0
	int loci0 = PotFFT.LocalI0();
	
	std::vector< std::complex<double> >
		X(PotFFT.LocalDim());

	int dim1=10000, dim2=1, norb =4; 
	for( int n=0 ; n< locDim0*dim1*dim2 ; n++)
	for( int o=0 ; o< norb ; o++)
		X[ n*norb + o ] = cos ( sqrt(o) *(n+loci0*dim1*dim2)*M_PI/3.0 );


	//Start random number variables
	gsl_rng_env_setup();
	
	const gsl_rng_type 
		*T = gsl_rng_default;

	gsl_rng
		*r = gsl_rng_alloc (T);

  int i, n = 10;

  for (i = 0; i < n; i++) 
    {
      double u = gsl_rng_uniform (r);
      printf ("%.5f\n", u);
    }

  gsl_rng_free (r);


  double packet_size;
  long seed;

  gsl_rng *rng;  // random number generator
  rng = gsl_rng_alloc (gsl_rng_default);     // pick random number generator
  seed = ext_seed* getpid()*( world_rank + 2321);    
  gsl_rng_set (rng, seed);                  // set seed

  packet_size = gsl_rng_uniform (1500); // get a random number from
                                            // the exponential distribution
  gsl_rng_free (rng);                       // dealloc the rng   

  return 0;

//	if ( world_rank == 0)
//	{
//		std::cout<<" Input "<<std::endl;
//		int o=1 ;
//		for( int n=0 ; n< 10 ; n++)
//			std::cout<< (n+loci0*dim1*dim2)*norb + o  <<" "<<X[ n*norb + o ] <<std::endl;
//		std::cout<<std::endl<<std::endl;
//	}
	int Mom=1; Mom=1;
 
 	PotFFT.DirectFourierTransform ( &X[0] );
	PotFFT.InverseFourierTransform( &X[0] );

//	if ( world_rank == 0)
//	{
//		std::cout<<" Output "<<std::endl;
//		int o=1 ;
//		for( int n=0 ; n< 10 ; n++)
//			std::cout<< (n+loci0*dim1*dim2)*norb + o  <<" "<<X[ n*norb + o ] <<std::endl;
//		std::cout<<std::endl<<std::endl;
//	}


	MPI_Finalize();
return 0;}
