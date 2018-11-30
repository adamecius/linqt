#ifndef PROJ_FFTW3_MPI_HPP
#define PROJ_FFTW3_MPI_HPP

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


class BZProjection
{

	public: 
	//Constructor
	BZProjection( ):rank(0)
	{
		in_= NULL;
		out_= NULL;
		filterIdx_=NULL;
	};
			
	//Methods
	void CleanUp()
	{
		fftw_destroy_plan(planForw_);
		fftw_destroy_plan(planBack_);
		
		if(in_!= NULL) fftw_free(in_); 
		if(out_!= NULL) fftw_free(out_);
		if(filterIdx_!=NULL) free(filterIdx_);
		fftw_cleanup();
	}

	void StartProjection(MPI_Comm& comm, std::string _projLabel )
	{
		// Get the mpi information that is going to be used in the class
		int world_size;
		MPI_Comm_size(comm, &world_size);
		MPI_Comm_rank(comm, &world_rank);
		MPI_Get_processor_name(processor_name, &name_len);
		
		//Tries to open the file projFileName_  
		projLabel_=_projLabel;
		std::string projFileName_="operators/"+projLabel_+".KPROJ";
		std::ifstream projFile;
		projFile.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		try  {	projFile.open ( projFileName_.c_str() );	}
		catch (std::ifstream::failure e)
		{
			projFile.close();
			if ( world_rank == 0)
			std::cerr 	<<"Exception:"<<e.what()
						<<"\n while opening config file: "
						<< projFileName_<<"."<<std::endl;
		}
		// If successful, the read the dimensions of the system
		// and the number of sites
		ptrdiff_t orbNum, maxSpin, orbPerCell, numCells;
		projFile>>n[0]>>n[1]>>n[2]>> orbNum>> maxSpin;	
		numCells 	= n[0]*n[1]*n[2];
		orbPerCell 	= orbNum*maxSpin;
		dim_		= numCells*orbPerCell;
		if( SpDim() == 0)
		{ 
			SetSpDim(3);	
			std::cout 	<<"The spatial dimension "
						<<"was not set, so we used "<<SpDim()<<"."<<std::endl;
		}
		else
			std::cout 	<<"The spatial dimension  was set to "<<SpDim()<<"."<<std::endl;
			
			
		// Automatically determine the  local data size of memory
		// depending on the dimensionality
		ptrdiff_t rank = SpDim();
		local_alloc =  fftw_mpi_local_size_many(rank, n, 
										orbPerCell,FFTW_MPI_DEFAULT_BLOCK,
										comm, &local_dim0, &local_i0);
		const ptrdiff_t //local dimension of the subspace per process 
		SetLocalDim(local_dim0*n[1]*n[2]*orbPerCell);

		// Allocate the space for the input and output vectors										
		in_  = (fftw_complex*) fftw_alloc_complex(local_alloc);
		out_ = (fftw_complex*) fftw_alloc_complex(local_alloc);
		filterIdx_ = (double*) malloc( sizeof(double) * LocalDim() );


		
//		std::string planLabel=(projLabel_+".MPI_FFTF");
//{		
		//Choose the best plans for the particular system
		//		bool foundPlan = fftw_import_wisdom_from_filename(planLabel.c_str()) ;
		// export	
//		{
//			int rank;
//			fftw_mpi_init();
//			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//			if (rank == 0) fftw_import_wisdom_from_filename("mywisdom");
//			fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD);
//		}

//		import
//		(Note that we must call fftw_mpi_init before importing any wisdom that might contain MPI plans.) Similarly, a typical code snippet to export wisdom from all processes to a file is:
//		{
//			int rank;
//
//		fftw_mpi_gather_wisdom(MPI_COMM_WORLD);
//			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//			if (rank == 0) fftw_export_wisdom_to_filename("mywisdom");
//		}
//		
//}		
		planForw_ = fftw_mpi_plan_many_dft(rank, n, orbPerCell, 
							FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, 
							in_,out_, comm,FFTW_FORWARD,FFTW_PATIENT);

		planBack_ = fftw_mpi_plan_many_dft(rank, n, orbPerCell, 
							FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, 
							in_,out_, comm,FFTW_BACKWARD,FFTW_PATIENT);

		//if(!foundPlan)
		//	fftw_export_wisdom_to_filename(planLabel.c_str());
//		projFile.exceptions ( std::ifstream::badbit );
//		while ( !projFile.eof() )
//		{
//			int n0;
//			double Idx;
//			projFile>>n0>>Idx;
//			for( int orb= 0 ; orb< orbPerCell&& !projFile.eof() ; orb ++ )
//			filterIdx_[n0*orbPerCell+orb]=Idx;
//		}
//		projFile.close();						

		//for testing purposes
//						for( int n0=0; n0< numCells; n0++)
//						for( int orb= 0 ; orb< orbPerCell ; orb ++ )
//							filterIdx_[n0*orbPerCell+orb]=1.0; 
												
//				};

	};	


/*

	void DirectFourierTransform( int size, std::complex<double>*   X	)
	{
		memcpy( &in_[0] , &X[0]   , LocalDim()*sizeof( fftw_complex ) );//Pass the vector to init
		fftw_execute( planForw_ ,&in_[0], &out_[0]); //Transfor the vector from real to momentum space
		memcpy( &X[0]   , &out_[0], Dim()*sizeof( fftw_complex ) );//Pass Final result to X
	}

	void InverseFourierTransform( int size,std::complex<double>*   X)
	{
		memcpy( &in_[0] , &X[0]   , LocalDim()*sizeof( fftw_complex ) );//Pass the vector to init
		fftw_execute( planForw_ ,&in_[0], &out_[0]); //Transfor the vector from real to momentum space
		memcpy( &X[0]   , &out_[0], Dim()*sizeof( fftw_complex ) );//Pass Final result to X
	}
			

	void ProjectVector( int size,std::complex<double>*   X)
	{
		if (size != Dim() )
			std::cerr<<"The dimension of the input vetor in "
					 <<"ProjectVector method of Projection class does "
					 <<"not mathc internal dimension"<<std::endl;
	
		memcpy( &in_[0], &X[0]   , Dim()*sizeof( fftw_complex ) );//Pass the vector to init		
		fftw_execute( planForw_ ); //Transfor the vector from real to momentum space
		for( int i=0;i < Dim(); i++) //Apply projection
		{	
			(out_[i])[0] = (out_[i])[0]*filterIdx_[i];
			(out_[i])[1] = (out_[i])[1]*filterIdx_[i];
		}		
		fftw_execute( planBack_ ); //Transfor the vector from momentum to real space
		memcpy( &X[0], &in_[0] , Dim()*sizeof( fftw_complex ) );
 
		//Pass Final result to X and normalize
		double numCells = n[0]*n[1]*n[2];
		for(int i=0; i< Dim() ;i++)
			X[i]/=numCells;


	}
	
*/
	//Returns the total dimension of the system
	ptrdiff_t LocalDim() const { return local_dim_; }
	ptrdiff_t Dim() const { return dim_; }
	ptrdiff_t SpDim() const { return rank; }
	void SetSpDim(const ptrdiff_t _rank) { rank = _rank; }
	void SetLocalDim(const ptrdiff_t _local_dim)  { local_dim_=_local_dim; }


	private: 
	ptrdiff_t rank ;
	ptrdiff_t dim_,local_dim_;
	ptrdiff_t  n[3];
	fftw_complex *in_,*out_;
	double* filterIdx_;
	fftw_plan planForw_, planBack_ ;
	std::string projLabel_;

	//mpi information for the class
	int world_size, world_rank, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
	ptrdiff_t local_alloc, local_dim0, local_i0;

};
#endif
