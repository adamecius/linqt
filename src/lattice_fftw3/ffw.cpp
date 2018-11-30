

#include "lattice_fftw3_mpi.hpp"

//The advanced-interface fftw_mpi_plan_many_dft additionally allows you 
//to specify the block sizes for the first dimension (block) of the 
//n0 × n1 × n2 × … × nd-1 input data and the first dimension (tblock) 
//of the n1 × n0 × n2 ×…× nd-1 transposed data (at intermediate steps of 
//the transform, and for the output if
// FFTW_TRANSPOSED_OUT is specified in flags). 
//These must be the same block sizes as were passed to the corresponding
// ‘local_size’ function; you can pass FFTW_MPI_DEFAULT_BLOCK to
// use FFTW’s default block size as in the basic interface. Also, the 
//howmany parameter specifies that the transform is of contiguous 
//howmany-tuples rather than individual complex numbers; this 
//corresponds to the same parameter in the serial advanced interface 
//(see Advanced Complex DFTs) with stride = howmany and dist = 1.

void LatticeFFTW::DirectFourierTransform( qt::complex*   X	)
{
	memcpy( &in_[0] , &X[0]   , LocalDim()*sizeof( fftw_complex ) );//Pass the vector to init
	fftw_mpi_execute_dft( planForw_ ,&in_[0], &out_[0]); //Transfor the vector from real to momentum space
	memcpy( &X[0]   , &out_[0], LocalDim()*sizeof( fftw_complex ) );//Pass Final result to X

}

void LatticeFFTW::InverseFourierTransform( qt::complex*   X)
{
	memcpy( &in_[0] , &X[0]   , LocalDim()*sizeof( fftw_complex ) );//Pass the vector to init
	fftw_mpi_execute_dft( planBack_ ,&in_[0], &out_[0]); //Transfor the vector from real to momentum space
	memcpy( &X[0]   , &out_[0], LocalDim()*sizeof( fftw_complex ) );//Pass Final result to X

//Pass Final result to X and normalize

}
	
