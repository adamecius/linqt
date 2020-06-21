// Used for OPENMP functions
#include "chebyshev_solver.hpp"


using namespace chebyshev;


std::array<double,2> utility::SpectralBounds( SparseMatrixType& HAM)
{
	double highest=0, lowest=0;
	std::ifstream bounds_file( "BOUNDS" );
	if( bounds_file.is_open() )
	{
		bounds_file>>lowest>>highest;
		bounds_file.close();
	}
	else
	{
		lowest = -100;
		highest = 100;
	}
	return { lowest, highest};
};


int chebyshev::ComputeMomTable( chebyshev::Vectors &chevVecL, chebyshev::Vectors& chevVecR ,  chebyshev::Moments2D& sub)
{
	const size_t maxMR = sub.HighestMomentNumber(0);
	const size_t maxML = sub.HighestMomentNumber(1);
		
	const int nthreads = mkl_get_max_threads();
//	mkl_set_num_threads_local(1); 
	for( auto m0 = 0; m0 < maxML; m0++)
	{
		//#pragma omp parallel for default(none) shared(chevVecL,chevVecR,m0,sub)
		for( auto m1 = 0; m1 < maxMR; m1++)
			sub(m0,m1) = linalg::vdot( chevVecL.Vector(m0) , chevVecR.Vector(m1) );
	}
//	mkl_set_num_threads_local(nthreads); 

return 0;
};


int chebyshev::CorrelationExpansionMoments(	const vector_t& PhiR, const vector_t& PhiL,
											SparseMatrixType &OPL,
											SparseMatrixType &OPR,  
											chebyshev::Vectors &chevVecL,
											chebyshev::Vectors &chevVecR,
											chebyshev::Moments2D &chebMoms
											)
{
	const size_t numRVecs = chevVecR.NumberOfVectors();
	const size_t numLVecs = chevVecL.NumberOfVectors();
	const size_t NumMomsR = chevVecR.HighestMomentNumber();
	const size_t NumMomsL = chevVecL.HighestMomentNumber();
	const size_t momvecSize = (size_t)( numRVecs*numLVecs );

	auto start = std::chrono::high_resolution_clock::now();
	chebyshev::Moments2D sub(numLVecs,numRVecs);

	
	chevVecL.SetInitVectors( OPL, PhiL );
	for(int  mL = 0 ; mL <  NumMomsL ; mL+=numLVecs )
	{
		chevVecL.IterateAll();

		chevVecR.SetInitVectors( PhiR );
		for(int  mR = 0 ; mR <  NumMomsR ; mR+=numRVecs)
		{
			chevVecR.IterateAll();
			chevVecR.Multiply( OPR );
			chebyshev::ComputeMomTable(chevVecL,chevVecR, sub );		
			chebMoms.AddSubMatrix(sub,mL,mR);
		}
	}	   
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Calculation of the moments in sequential solver took: " << elapsed.count() << " seconds\n\n";

	return 0;
};



int chebyshev::CorrelationExpansionMoments( SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::Moments2D &chebMoms, qstates::generator& gen )
{
	//Set Batch Behavior
	int batchSize;
    if( !chebyshev::GetBatchSize(batchSize)  )
    {
		batchSize = chebMoms.HighestMomentNumber(0);
		std::cout<<"Using default BATCH_SIZE = "<<batchSize<<std::endl;
	}
	else
		std::cout<<"Using user's BATCH_SIZE = "<<batchSize<<std::endl;


	int 
	or_mom0=chebMoms.HighestMomentNumber(0),
	or_mom1=chebMoms.HighestMomentNumber(1);
	if( or_mom0%batchSize != 0 ||
		or_mom1%batchSize != 0 )
	{
		std::cout<<"\nWARNING: This implementation need to the batchsize to be a multiple of ";
		std::cout<<or_mom0<<" and "<<or_mom1<<std::endl;
		std::cout<<"The moments requiremenets will be increased to ";
		
		int mom0= ( (or_mom0+batchSize-1)/batchSize )*batchSize;
		int mom1= ( (or_mom1+batchSize-1)/batchSize )*batchSize;
		std::cout<<mom0<<"x"<<mom1<<std::endl;
		std::cout<<" which is suboptimal "<<std::endl<<std::endl;
		chebMoms = chebyshev::Moments2D( chebMoms, mom0,mom1);
	}

	
	const int DIM  = chebMoms.SystemSize(); 
	
	chebyshev::Vectors 
	chevVecL( chebMoms ), chevVecR( chebMoms );

	//Allocate the memory
	
	chevVecL.SetNumberOfVectors( batchSize );
	chevVecR.SetNumberOfVectors( batchSize );
	printf("Chebyshev::CorrelationExpansionMoments will used %f GB\n", chevVecL.MemoryConsumptionInGB() + chevVecR.MemoryConsumptionInGB() );

	//This operation is memory intensive
	std::cout<<"Initializing chevVecL"<<std::endl;
	chevVecL.CreateVectorSet( );
	std::cout<<"Initialize chevVecR"<<std::endl;
	chevVecR.CreateVectorSet( );

	gen.SystemSize(DIM);
	while( gen.getQuantumState() )
	{
		std::cout<<"Computing with ID: "<<gen.count<<" states" <<std::endl;

		//SELECT RUNNING TYPE
		chebyshev::CorrelationExpansionMoments(gen.State(),gen.State(), OPL, OPR, chevVecL,chevVecR, chebMoms);
	}

	//Fix the scaling of the moments
    const int NumMomsL = chebMoms.HighestMomentNumber(0);
    const int NumMomsR = chebMoms.HighestMomentNumber(1);
	for (int mL = 0 ; mL < NumMomsL; mL++)				  
	for (int mR = mL; mR < NumMomsR; mR++)
	{
		double scal=4.0/gen.NumberOfStates();
		if( mL==0) scal*=0.5;
		if( mR==0) scal*=0.5;

		const value_t tmp = scal*( chebMoms(mL,mR) + std::conj(chebMoms(mR,mL)) )	/2.0;
		chebMoms(mL,mR)= tmp;
		chebMoms(mR,mL)= std::conj(tmp);
	}
	
	
	if( NumMomsL != or_mom0 or NumMomsR != or_mom1 )
	{
		std::cout<<"Changin from "<<NumMomsL<<"x"<<NumMomsR<<" to "<<or_mom0<<"x"<<or_mom1<<std::endl;
		chebMoms.MomentNumber(or_mom0,or_mom1);
	}
	 
return 0;
};



int chebyshev::SpectralMoments( SparseMatrixType &OP,  chebyshev::Moments1D &chebMoms, qstates::generator& gen )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();

	gen.SystemSize(Dim);
	while( gen.getQuantumState() )
	{		
		auto Phi = gen.State();
		//Set the evolved vector as initial vector of the chebyshev iterations
		if (OP.isIdentity() )
			chebMoms.SetInitVectors( Phi );
		else
			chebMoms.SetInitVectors( OP,Phi );
			
		for(int m = 0 ; m < NumMoms ; m++ )
		{
			double scal=2.0/gen.NumberOfStates();
			if( m==0) scal*=0.5;
			chebMoms(m) += scal*linalg::vdot( Phi, chebMoms.Chebyshev0() ) ;
			chebMoms.Iterate();
		}
	}
	return 0;
};


int chebyshev::TimeDependentCorrelations(SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::MomentsTD &chebMoms, qstates::generator& gen  )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	const auto NumTimes= chebMoms.MaxTimeStep();
	
	//Initialize the Random Phase vector used for the Trace approximation
	gen.SystemSize(Dim);	
	while( gen.getQuantumState() )
	{
		chebMoms.ResetTime();

		auto PhiR =gen.State();
		auto PhiL =gen.State();
		 
		//Multiply right operator its operator
		OPL.Multiply(PhiR,PhiL); //Defines <Phi| OPL 
		
		//Evolve state vector from t=0 to Tmax
		while ( chebMoms.CurrentTimeStep() !=  chebMoms.MaxTimeStep()  )
		{
			const auto n = chebMoms.CurrentTimeStep();

			//Set the evolved vector as initial vector of the chebyshev iterations
			chebMoms.SetInitVectors( OPR , PhiR );

			for(int m = 0 ; m < NumMoms ; m++ )
			{
				double scal=2.0/gen.NumberOfStates();
				if( m==0) scal*=0.5;
				chebMoms(m,n) += scal*linalg::vdot( PhiL, chebMoms.Chebyshev0() ) ;
				chebMoms.Iterate();
			}
			
			chebMoms.IncreaseTimeStep();
			//evolve PhiL ---> PhiLt , PhiR ---> PhiRt 
			chebMoms.Evolve(PhiL) ;
			chebMoms.Evolve(PhiR) ;
		}
	
	}
	
	return 0;
};
