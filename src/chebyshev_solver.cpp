// Used for OPENMP functions
#include "chebyshev_solver.hpp"


using namespace chebyshev;



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



int parallel::CorrelationExpansionMoments( const vector_t& PhiR, const vector_t& PhiL,
											SparseMatrixType &HAM,
											SparseMatrixType &OPL,
											SparseMatrixType &OPR,  chebyshev::Moments2D &chebMoms)
{
    const int DIM = HAM.rank();
	const int NumMomsR= chebMoms.HighestMomentNumber(0);
	const int NumMomsL= chebMoms.HighestMomentNumber(1);

	auto start = std::chrono::high_resolution_clock::now();

	chebyshev::Vectors chevVecL( chebMoms,0 );
	chebyshev::Vectors chevVecR( chebMoms,1 );
	int batchSize = NumMomsR;
	chebyshev::Moments::vector_t momvec( batchSize*batchSize );


	printf("hebyshev::parallel::CorrelationExpansionMoments will used %f GB\n",chevVecL.MemoryConsumptionInGB()+chevVecR.MemoryConsumptionInGB());

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
	assert( chebVL.SystemSize() == chebVR.SystemSize() && chebVL.HighestMomentNumber() == chebVR.HighestMomentNumber() );
	const int dim    = chebVL.SystemSize();
	const int numMom = chebVL.HighestMomentNumber() ;
	assert( output.size() == numMom*numMom );

	linalg::batch_vdot(dim,numMom,&chebVL(0),&chebVR(0),&output[0]);
return 0;
};


int chebyshev::CorrelationExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::Moments2D &chebMoms, StateType type )
{

	//Allocate the memory
    string NUM_THREADS ="Default"; 	if(getenv("OMP_NUM_THREADS")) NUM_THREADS = getenv("OMP_NUM_THREADS");
	int batchSize=3;
    bool use_sequential = !( chebyshev::GetBatchSize(batchSize) );


	std::cout<<"The Correlation calculation will run on "<<NUM_THREADS<< " threads"<<std::endl;
	chebMoms.Rescale2ChebyshevDomain(HAM);

    const int DIM = HAM.rank(); 

	//allocate the memory for the input vector, and the iteration vector
	std::vector< std::complex<double> >  Phi(DIM); 	//States Vectors
	for (int i = 0; i < numStates; i++)
	{
		//SELECT STATE TYPE
		std::cout<<"Computing state: "<<i+1<<" of "<<numStates<<std::endl;

		for(int j = 0; j < DIM ; j++)
		switch (type)
		{
			case RANDOM_STATE:
			Phi[j] = exp(value_t(0, 2.0*M_PI*(double)rand() / (double)RAND_MAX ) )/sqrt(DIM);
			break;
			case LOCAL_STATE:
			Phi[j] = ( (j==i) ? 1.0 : 0.0 ) ;
			break;
			default:
			std::cerr<<" The state state is not identify, aborting running"<<std::endl;
			std::exit(-1);
		}


		//SELECT RUNNING TYPE
		if( use_sequential )
		{
			std::cout<<"USING SEQUENTIAL IMPLEMENTATION"<<std::endl;
			chebyshev::sequential::CorrelationExpansionMoments(	Phi,Phi, HAM, OPL, OPR, chebMoms);
		}
		else
		{
			std::cout<<"LARGE MEMORY CALCULATION"<<std::endl;
			chebyshev::parallel::CorrelationExpansionMoments( Phi,Phi, HAM, OPL, OPR, chebMoms);
		}
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


