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
		std::cout<<"File BOUNDS not found, computing spectral boundas automatically (EXPERIMENTAL)"<<std::endl;
		const int NUMIT = 100;
		const int DIM = HAM.rank(); 
		double absmax=0;
		vector_t PhiL(DIM), PhiR(PhiL);
		for( int num_extrema = 0 ; num_extrema < 2  ; num_extrema++)
		{
			qstates::FillWithRandomPhase(PhiL);
			for(int it =0; it  < NUMIT ; it++)
			{
				HAM.Multiply(PhiL, PhiR);
				absmax = linalg::nrm2(PhiR); ////Get <Phi|H^2|Phi>
				linalg::scal(1.0/absmax, PhiR);				
				
				//Break is goal achieved
				if( std::fabs(highest-absmax) < 1e-2 )
				{
					HAM.Multiply(PhiR, PhiL);
					highest = linalg::vdot(PhiL,PhiR).real(); ////Get <Phi|H^2|Phi>	
					break;
				}
				//else continue
				linalg::copy(PhiR,PhiL);
				highest = absmax;
			}
			if( lowest == 0)
			{
				HAM.Rescale(1.0,-highest);					
				lowest = highest;
				highest = 0;
			}
		}
		if( lowest > highest )
		{
			absmax = highest;
			highest= lowest ;
			lowest = absmax;
		}
	}
	std::cout<<lowest<<" "<<highest<<" "<<(highest-lowest)<<" "<<(lowest+highest)/2<<std::endl;
	return { lowest, highest};

};

int sequential::KuboGreenwoodChebMomConvergence( const double E0,
												 const double eta,
												 SparseMatrixType &HAM,
												 SparseMatrixType &OPL,
												 SparseMatrixType &OPR,
												 chebyshev::Moments1D&  chebMoms)
{
	const int numChebVec  = 1 , DIM = chebMoms.SystemSize(); 
	const int MOM  = chebMoms.HighestMomentNumber();
	const int deltaMOM  = chebMoms.JacksonKernelMomCutOff(eta);
	std::cout<<"deltaMom "<<deltaMOM<<std::endl;
	//Initialize chebyshev vectors
	chebyshev::Vectors chebVecs(numChebVec,DIM);

	//Normalize the hamiltonian, energies, and get geometric factors
	const double geo_fact =  DIM/chebMoms.HalfWidth()/chebMoms.HalfWidth();
	const double energ    =  chebMoms.Rescale2ChebyshevDomain(E0);
	chebMoms.Rescale2ChebyshevDomain(HAM);
	
	//Create a random phase vector
	vector_t PhiR(DIM),PhiL(DIM);
	qstates::FillWithRandomPhase(PhiR); 

	//COMPUTE VECTOR <Phi OPR delta(H-E)| 
	chebVecs.SetInitVectors( HAM, OPR, PhiR ); //<j0|= <Phi|OPR
	linalg::scal( 0.0, PhiL );
	for( int m=0; m < deltaMOM; m++)
	{
		auto chebCL = chebMoms.JacksonKernel(m,deltaMOM)*delta_chebF(energ,m); if(m==0) chebCL*=0.5; 
		linalg::axpy( chebCL, chebVecs.Chebyshev0(),  PhiL);
		chebVecs.Iterate( HAM );
	}
	std::cout<<"Finished adding broadening of "<<eta<<std::endl;
	//COMPUTE VECTOR <Phi OPR delta(H-E)  |OPL G | Phi>  
	chebVecs.SetInitVectors( HAM, OPL, PhiL ); //<j0|= <Phi OPR delta(H-E)  |OPL
	linalg::scal( 0.0, PhiL );
	for( int m=0; m < MOM; m++)
	{
		auto chebCR = greenR_chebF(energ,m); if(m==0) chebCR*=0.5; 
		linalg::axpy( chebCR, chebVecs.Chebyshev0(),  PhiL);
		chebMoms(m) = linalg::vdot( PhiL, PhiR ).imag()*geo_fact; //This actually gives <JR|JL>*
		chebVecs.Iterate( HAM );

		std::cout<<m<<" "<<chebMoms(m).real()<<" "<<std::endl;
	}

	return 0;
};							  


int sequential::DensityExpansionMoments(vector_t& PhiL,vector_t& PhiR,
							SparseMatrixType &HAM,
							SparseMatrixType &OP,
							chebyshev::Moments1D &chebMoms)
{ 
    const int DIM  = chebMoms.SystemSize(); 
	const double A = chebMoms.ScaleFactor();
    const double B = chebMoms.ShiftFactor();
    const int MOM  = chebMoms.HighestMomentNumber();

	vector_t J0(DIM), J1(DIM), JOPL(DIM);
	
	linalg::copy(PhiL,J0);
	HAM.Multiply(J0,J1);
	for(int m=0; m < MOM ; m++)
	{
		OP.Multiply(J0, JOPL);
		chebMoms(m) += linalg::vdot(JOPL,PhiR); //This actually gives <JR|JL>*
		HAM.Multiply(2.0,J1,-1.0,J0);
		J0.swap(J1);
	}
	return 0;
};

int sequential::ComputeMomTable( chebyshev::Vectors &chebVL, chebyshev::Vectors& chebVR ,  vector_t& output)
{
	std::cout<<"NOT IMPLEMENTED sequential::ComputeMomTable"<<std::endl;
return 0;
};


int sequential::CorrelationExpansionMoments( const vector_t& PhiL,const vector_t& PhiR,
											  SparseMatrixType &HAM,
											  SparseMatrixType &OPL,
											  SparseMatrixType &OPR,
											  chebyshev::Moments2D &chebMoms)
{
	const int DIM  = chebMoms.SystemSize(); 
	const double A = chebMoms.ScaleFactor();
    const double B = chebMoms.ShiftFactor();
    const int MOML = chebMoms.HighestMomentNumber(0);
    const int MOMR = chebMoms.HighestMomentNumber(1);
    
    vector_t JR0(DIM),JR1(DIM),JL0(DIM),JL1(DIM),JOL(DIM);

    //Start the chebyshev expansion of the correlations
	auto start = std::chrono::high_resolution_clock::now();
	OPR.Multiply(PhiR, JR0);
    HAM.Multiply(JR0 , JR1);
    for (int m1 = 0; m1 < MOMR; m1++ )
    {
		linalg::copy(PhiL,JL0);
		HAM.Multiply( JL0, JL1 );
		for (int m0= 0; m0 < MOML; m0++)
		{
			OPL.Multiply(JL0, JOL );
			chebMoms(m0,m1) += linalg::vdot(JOL,JR0); //This actually gives <JR|JL>*
			HAM.Multiply(2.0, JL1, -1.0, JL0);
			JL0.swap(JL1);
		};
		HAM.Multiply(2.0, JR1, -1.0, JR0 );
		JR0.swap(JR1);
	}
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Calculation of the moments in sequential solver took: " << elapsed.count() << " seconds\n";
	
	return 0;
};



int parallel::CorrelationExpansionMoments(	const int batchSize,
											const vector_t& PhiR, const vector_t& PhiL,
											SparseMatrixType &HAM,
											SparseMatrixType &OPL,
											SparseMatrixType &OPR,  
											chebyshev::Vectors &chevVecL,
											chebyshev::Vectors &chevVecR,
											chebyshev::Moments2D &chebMoms
											)
{
	const size_t DIM = HAM.rank();
	const size_t NumMomsR = chevVecR.HighestMomentNumber();
	const size_t NumMomsL = chevVecL.HighestMomentNumber();
	const size_t momvecSize = (size_t)( (long unsigned int)batchSize*(long unsigned int)batchSize );


	auto start = std::chrono::high_resolution_clock::now();
	
	std::cout<<"Initialize sparse for moment matrix"<<std::endl;
	chebyshev::Moments::vector_t momvec( momvecSize );


	for(int  mR = 0 ; mR <  NumMomsR ; mR+=batchSize)
	{
		chevVecL.SetInitVectors( HAM, OPL, PhiL );
		chevVecL.IterateAll( HAM );
		for(int  mL = 0 ; mL <  NumMomsL ; mL+=batchSize)
		{
			chevVecR.SetInitVectors( HAM, PhiR );
			chevVecR.IterateAll( HAM );
			chevVecR.Multiply( OPR );
			parallel::ComputeMomTable(chevVecL,chevVecR, momvec );		
			linalg::axpy(momvec.size(), 1.0 , &momvec[0], &chebMoms(mR,mL) );
		}
	}	   
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Calculation of the moments in sequential solver took: " << elapsed.count() << " seconds\n";
	
	return 0;
};


int parallel::ComputeMomTable( chebyshev::Vectors &chebVL, chebyshev::Vectors& chebVR ,  vector_t& output)
{
	const auto dim   = chebVL.SystemSize();
	const auto maxML = chebVL.HighestMomentNumber();
	const auto maxMR = chebVR.HighestMomentNumber();
	assert( chebVL.SystemSize() == chebVR.SystemSize() && chebVL.HighestMomentNumber() == chebVR.HighestMomentNumber() );
	assert( output.size() == maxML*maxMR );
	const int nthreads = mkl_get_max_threads();

	mkl_set_num_threads_local(1); 
	for( auto m0 = 0; m0 < maxML; m0++)
	{
		#pragma omp parallel for default(none) shared(chebVL,chebVR,m0,output)
		for( auto m1 = 0; m1 < maxMR; m1++)
			output[ m0*maxMR + m1 ] = linalg::vdot( chebVL.Vector(m0) , chebVR.Vector(m1) );	
	}
	mkl_set_num_threads_local(nthreads); 

//	linalg::batch_vdot(dim,numMom,&chebVL(0),&chebVR(0),&output[0]);
return 0;
};


int chebyshev::CorrelationExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::Moments2D &chebMoms, StateType type )
{
    int kpm_seed = time(0); 	if(getenv("KPM_SEED")) kpm_seed = std::stoi(string(getenv("KPM_SEED")));
	srand(kpm_seed);
	std::cout<<"Using seed "<<kpm_seed<<std::endl;
	//Allocate the memory
    string NUM_THREADS ="Default"; 	if(getenv("OMP_NUM_THREADS")) NUM_THREADS = getenv("OMP_NUM_THREADS");
	int batchSize=3;
    bool use_sequential = !( chebyshev::GetBatchSize(batchSize) );


	std::cout<<"The Correlation calculation will run on "<<NUM_THREADS<< " threads"<<std::endl;
	chebMoms.Rescale2ChebyshevDomain(HAM);
    const int DIM = HAM.rank(); 


	chebyshev::Vectors chevVecL,chevVecR;
	if( use_sequential )
		std::cout<<"USING SEQUENTIAL IMPLEMENTATION"<<std::endl;
	else
	{
		std::cout<<"USING LARGE MEMORY IMPLEMENTATION with BATCH_SIZE: "<<batchSize<<std::endl;
		printf("Chebyshev::parallel::CorrelationExpansionMoments will used %f GB\n", chevVecL.MemoryConsumptionInGB() + chevVecR.MemoryConsumptionInGB() );

		std::cout<<"Initialize chevVecL"<<std::endl;
		chevVecL=chebyshev::Vectors( chebMoms,0 );
		std::cout<<"Initialize chevVecR"<<std::endl;
		chevVecR=chebyshev::Vectors( chebMoms,1 );
	}

	//allocate the memory for the input vector, and the iteration vector
	std::vector< std::complex<double> >  Phi(DIM); 	//States Vectors
	for (int i = 0; i < numStates; i++)
	{

		for(int j = 0; j < DIM ; j++)
		switch (type)
		{
			case RANDOM_STATE:
			Phi[j] = exp(value_t(0, 2.0*M_PI*(double)rand() / (double)RAND_MAX ) )/sqrt(DIM);
			if( j <10)
				std::cout<<j<<" "<<Phi[j]<<std::endl;

			break;
			case LOCAL_STATE:
 //			i=(numStates-1); 
			Phi[j] = ( (j==i) ? 1.0 : 0.0 ) ;
			break;
			default:
			std::cerr<<" The state state is not identify, aborting running"<<std::endl;
			std::exit(-1);
		}
//		numStates= 1;
		//SELECT STATE TYPE
		std::cout<<"Computing with ID: "<<i+1<<" of "<<numStates<<" states" <<std::endl;

		//SELECT RUNNING TYPE
		if( use_sequential )
			chebyshev::sequential::CorrelationExpansionMoments(	Phi,Phi, HAM, OPL, OPR, chebMoms);
		else
			chebyshev::parallel::CorrelationExpansionMoments(batchSize, Phi,Phi, HAM, OPL, OPR, chevVecL,chevVecR, chebMoms);
	}

	//Fix the scaling of the moments
    const int NumMomsL = chebMoms.HighestMomentNumber(0);
    const int NumMomsR = chebMoms.HighestMomentNumber(1);
	for (int mL = 0 ; mL < NumMomsL; mL++)				  
	for (int mR = mL; mR < NumMomsR; mR++)
	{
		double scal=4.0/numStates;
		if( mL==0) scal*=0.5;
		if( mR==0) scal*=0.5;

		const value_t tmp = scal*( chebMoms(mL,mR) + std::conj(chebMoms(mR,mL)) )	/2.0;
		chebMoms(mL,mR)= tmp;
		chebMoms(mR,mL)= std::conj(tmp);
	}

return 0;
};


int chebyshev::SpectralMoments(int numStates, SparseMatrixType &OP,  chebyshev::Moments1D &chebMoms, StateType type )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	vector_t Phi(Dim);
	for(int s=0; s < numStates; s++)
	{
		std::cout<<"Computing state "<<s<<" for spectral moments"<<std::endl;
		//Initialize the Random Phase vector used for the Trace approximation
		qstates::FillWithRandomPhase(Phi);// Defines  |Psi>
		
		//Set the evolved vector as initial vector of the chebyshev iterations
		if (OP.isIdentity() )
			chebMoms.SetInitVectors( Phi );
		else
			chebMoms.SetInitVectors( OP,Phi );
			
		for(int m = 0 ; m < NumMoms ; m++ )
		{
			double scal=2.0/numStates;
			if( m==0) scal*=0.5;
			chebMoms(m) += scal*linalg::vdot( Phi, chebMoms.Chebyshev0() ) ;
			chebMoms.Iterate();
		}
	}
	
	return 0;
};


int chebyshev::TimeDependentCorrelations(int numStates, SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::MomentsTD &chebMoms, StateType type )
{
	const auto Dim = chebMoms.SystemSize();
	const auto NumMoms = chebMoms.HighestMomentNumber();
	const auto NumTimes= chebMoms.MaxTimeStep();
	
	//Initialize the Random Phase vector used for the Trace approximation
	vector_t PhiL(Dim), PhiR(Dim);
	
	for(int s=0; s < numStates; s++)
	{
		chebMoms.ResetTime();
		qstates::FillWithRandomPhase(PhiR);// Defines  |Psi>
			
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
				double scal=2.0/numStates;
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
