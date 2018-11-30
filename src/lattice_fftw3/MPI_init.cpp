

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

bool
LatticeFFTW::Init(MPI_Comm comm,std::vector<qt::dimension> _n,const qt::dimension orbPerCell, std::string projLabel )
{
	// Get the mpi information that is going to be used in the class
	MPI_Comm_size(comm, &world_size);
	MPI_Comm_rank(comm, &world_rank);
	MPI_Get_processor_name(processor_name, &name_len);

	//Set the system label and try to open the configuration file
	SetLabel( projLabel );
	// If successful, read the dimensions of the system
	// and the number of sites
	ptrdiff_t n[3], spatial=0;
	for(int i=0;i<3; i++ )
	{
		if(_n[i]== 0)
			return false;

		n[i] = _n[i];
		if( _n[i]!= 1 )
			spatial += 1;

	}
	// Define the internal and extenral degrees of freedom
	ptrdiff_t
	numCells 	= n[0]*n[1]*n[2],
	dim_		= numCells*orbPerCell;
	SetSpDim(spatial);
			
	// Automatically determine the  local data size of memory
	// depending on the dimensionality
	local_alloc =  	fftw_mpi_local_size_many(SpDim(), n,
					orbPerCell,FFTW_MPI_DEFAULT_BLOCK,
					comm, &local_dim0, &local_i0);
	//local dimension of the subspace per process 
	SetLocalDim( local_dim0*n[1]*n[2]*orbPerCell );

	if ( world_rank == 0)
	{
		std::cout	<<std::endl<<"FFTW_STATUS"<<std::endl
					<<"The dim0= "<<n[0]<<" was splitted into "
					<<"local_dim0="<<local_dim0<<std::endl;

		std::cout	<<"The number of elements to be used by FFT is: "
					<<local_alloc<<"\nWhile the local dimension is:"
					<<LocalDim()<<std::endl;

		std::cout	<<"The memory cost for a single MPI-FFT will be: "
					<<(double)local_alloc*(double)sizeof( fftw_complex )/double(1048576)
					<<" MB"<<std::endl<<std::endl;

	}
	
	// Allocate the space for the input and output vectors										
	in_  = (fftw_complex*) fftw_alloc_complex(local_alloc);
	out_ = (fftw_complex*) fftw_alloc_complex(local_alloc);

	//Look for a previous plan and if found broadcast it to all nodes
	{
		if (world_rank == 0)
		{
			std::cout<<"Attemping to read plan wisdow from : "<< FFTWPlanPath()<<std::endl;
			fftw_import_wisdom_from_filename( FFTWPlanPath() );
		}
			fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD);
	}

	// Use either the previous wiswom to create the plan
	// or test different alternatives through the flag PATIENT
	planForw_ = fftw_mpi_plan_many_dft(rank, n, orbPerCell, 
						FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, 
						in_,out_, comm,FFTW_FORWARD,FFTW_PATIENT);
	planBack_ = fftw_mpi_plan_many_dft(rank, n, orbPerCell, 
						FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, 
						in_,out_, comm,FFTW_BACKWARD,FFTW_PATIENT);

	// After finding the ideal plan for this particular 
	// system configuration, save it in FFTPlanPath
	{
		int rank;
		fftw_mpi_gather_wisdom(MPI_COMM_WORLD);
		if (world_rank == 0)
		{ 
			std::cout<<"Saving plan wisdow in : "<< FFTWPlanPath()<<std::endl;
			fftw_export_wisdom_to_filename(FFTWPlanPath());
		}
	}

return true;	};	
