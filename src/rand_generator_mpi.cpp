#include <mpi.h>
#include <cmath>
#include <cstdlib>     /* srand, rand */
#include <iostream>
#include <sstream>
#include <vector>
#include <complex>
#include <gsl/gsl_rng.h>
#include <sys/types.h>
#include <unistd.h>

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int ext_seed = 234231;
	
	int dim0 = 10000 , dim1=10000, dim2=1, norb =4 , locDim0=( dim0 +world_size -1)/world_size; 
	int loci0 = world_rank*locDim0;

	//Start random number variables
	gsl_rng_env_setup();
	gsl_rng *rng;  // random number generator
	rng = gsl_rng_alloc (gsl_rng_default);     // uses the default
	int seed = ext_seed* getpid()*( world_rank + 2321);    //create a seed that is different for process.
	gsl_rng_set (rng, seed);                  // set seed


	std::vector< std::complex<double> > X(locDim0);
		
	for( int n=0 ; n< locDim0*dim1*dim2 ; n++)
	for( int o=0 ; o< norb ; o++)
	{
		double theta_r= ( 2.0*gsl_rng_uniform(rng) - 1.0 )*M_PI;
		X[ n*norb + o ] = std::complex<double>( cos( theta_r ), sin( theta_r) );
	}
	
	std::stringstream ss;
	ss << seed;
	string str = ss.str();

	std::string outputRe ("random/MyRandomRe-withSeed"+str+".dat");
	for( int n=0 ; n< locDim0*dim1*dim2 ; n++)
	for( int o=0 ; o< norb ; o++)
	{
		outputRe<<X[ n*norb + o ].real();
	}		
	outputRe.close();

	std::string outputIm ("random/MyRandomIm-withSeed"+str+".dat");
	for( int n=0 ; n< locDim0*dim1*dim2 ; n++)
	for( int o=0 ; o< norb ; o++)
	{
		outputIm<<X[ n*norb + o ].imag();
	}		
	outputIm.close();

	gsl_rng_free (rng);                       // dealloc the rng   
	MPI_Finalize();
return 0;}
