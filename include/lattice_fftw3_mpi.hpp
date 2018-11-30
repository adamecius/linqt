#ifndef LATTICE_FFTW3_MPI_HPP
#define LATTICE_FFTW3_MPI_HPP

#include "types_definitions.hpp"
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>
#include <fftw3-mpi.h>



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


class LatticeFFTW
{


	public: 
	//Constructor
	LatticeFFTW( ):rank(0)
	{
		in_= NULL;
		out_= NULL;
		fftw_mpi_init();
	};
			
	//Methods
	void CleanUp()
	{
		fftw_destroy_plan(planForw_);
		fftw_destroy_plan(planBack_);
		
		if(in_!= NULL) fftw_free(in_); 
		if(out_!= NULL) fftw_free(out_);
		fftw_mpi_cleanup();
	}

	bool 
	Init(MPI_Comm comm,std::vector<qt::dimension> _n,
		const qt::dimension orbPerCell, std::string projLabel );

	void DirectFourierTransform( qt::complex*   X	);

	void InverseFourierTransform( qt::complex*   X) ;
	

	//PUBLIC GETTERS AND SETTERS
	public: 
	
	inline std::string 
	Label() const { return fftw_label_;}
	
	inline const char* 
	FFTWPath() const { return fftw_path_.c_str();}
	
	inline const char* 
	FFTWPlanPath() const { return fftw_plan_label_.c_str();}
	
	inline ptrdiff_t 
	LocalDim() const { return local_dim_; }
	
	inline ptrdiff_t 
	LocalDim0() const { return local_dim0; }
	
	inline ptrdiff_t 
	LocalI0() const { return local_i0; }
	
	inline ptrdiff_t 
	Dim() const { return dim_; }
	
	inline ptrdiff_t 
	SpDim() const { return rank; }
	
	inline void 
	SetSpDim(const ptrdiff_t _rank) { rank = _rank; }
	
	inline void 
	SetLocalDim(const ptrdiff_t _local_dim)  { local_dim_=_local_dim; }
	
	inline void 
	SetLabel(std::string _fftw_label)
	{ 
		fftw_label_=_fftw_label;
		fftw_path_ ="operators/"+fftw_label_+".LAT";
		fftw_plan_label_ = fftw_label_+".MPI_FFTWIS";
	}


	private: 
	ptrdiff_t rank ;
	ptrdiff_t dim_,local_dim_;
	ptrdiff_t  n[3];
	fftw_complex *in_,*out_;
	fftw_plan planForw_, planBack_ ;

	//String variable where paths and inputs and output files are defined
	std::string fftw_label_, fftw_path_, fftw_plan_label_;

	//mpi information for the class
	int world_size, world_rank, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
	ptrdiff_t local_alloc, local_dim0, local_i0;
	
};


//#include "lattice_fftw3/MPI_init.cpp"
//#include "lattice_fftw3/ffw.cpp"
#endif
