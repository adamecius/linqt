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

int sequential::NonEqConvergence( SparseMatrixType &HAM,
							  SparseMatrixType &OPL,
							  SparseMatrixType &OPR,
							  double eta, double E0)
{
	chebyshev::Vectors chebVL, chebVR;
	
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


int parallel::TempDensityExpansionMoments(const int batchSize,					 
					  const vector_t& Phi,
					  SparseMatrixType &HAM,
					  SparseMatrixType &OP,
					  chebyshev::Vectors &chebVecR, chebyshev::Vectors &chebVecL,
					  chebyshev::MomentsTD &chebMoms)
{
  assert(chebVecR.HighestMomentNumber() == chebVecL.HighestMomentNumber());
  const size_t DIM = HAM.rank();
  const size_t NumMoms = chebVecR.HighestMomentNumber();
  const size_t NumTimes = chebMoms.HighestTimeNumber();
  const size_t momvecSize = (size_t)( (long unsigned int)batchSize );

  auto start = std::chrono::high_resolution_clock::now();

  std::cout << "Initialize sparse for moment matrix" << std::endl;
  chebyshev::Moments::vector_t momvec( momvecSize );

  for (int m = 0; m < NumMoms; m+=batchSize)
    {
      chebVecR.SetInitVectors(HAM, Phi);
      chebVecR.IterateAll(HAM);
      chebVecR.Multiply(OP);

      chebVecL.SetInitVectors(HAM, Phi);
      for (int n = 0; n < NumTimes; n++)
	{
	  chebVecL.EvolveAll(HAM, chebMoms.TimeStep(), chebMoms.TimeCoeff());
	  parallel::ComputeTDMomTable(chebVecL, chebVecR, momvec);
	  linalg::axpy(momvec.size(), 1.0, &momvec[0], &chebMoms(m, n));
	}
    }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Calculation of the moments in parallel solver took: "
	    << elapsed.count() << " seconds. " << std::endl;

  return 0;
};


int parallel::ComputeTDMomTable(chebyshev::Vectors &chebVL, chebyshev::Vectors& chebVR, vector_t& output)
{
  const auto dim   = chebVL.SystemSize();
  const auto maxM = chebVL.HighestMomentNumber();
  assert( chebVL.SystemSize() == chebVR.SystemSize() && chebVL.HighestMomentNumber() == chebVR.HighestMomentNumber() );
  assert( output.size() == maxM );
  const int nthreads = mkl_get_max_threads();
  
  mkl_set_num_threads_local(1);
#pragma omp parallel for default(none) shared(chebVL,chebVR,output)
  for( auto m = 0; m < maxM; m++)
    {
      output[m] = linalg::vdot( chebVL.Vector(m) , chebVR.Vector(m));	
    }
  mkl_set_num_threads_local(nthreads);
  
  return 0;
};


int chebyshev::TempDensityExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OP, SparseMatrixType &PROJ,  chebyshev::MomentsTD &chebMoms, StateType type )
{
  int kpm_seed = time(0); 	if(getenv("KPM_SEED")) kpm_seed = std::stoi(string(getenv("KPM_SEED")));
  srand(kpm_seed);
  std::cout<<"Using seed "<<kpm_seed<<std::endl;
  //Allocate the memory
  string NUM_THREADS ="Default"; 	if(getenv("OMP_NUM_THREADS")) NUM_THREADS = getenv("OMP_NUM_THREADS");
  int batchSize=3;
  bool use_sequential = !( chebyshev::GetBatchSize(batchSize) );
  
  
  std::cout<<"The Temporal Density calculation will run on "<<NUM_THREADS<< " threads"<<std::endl;
  chebMoms.Rescale2ChebyshevDomain(HAM);
  const int DIM = HAM.rank(); 

  
  chebyshev::Vectors chebVecL,chebVecR;
  /*if( use_sequential )
    std::cout<<"USING SEQUENTIAL IMPLEMENTATION"<<std::endl;
    else
    {*/
  assert( !(use_sequential) ); 
  std::cout<<"USING LARGE MEMORY IMPLEMENTATION with BATCH_SIZE"<<batchSize<<std::endl;
  printf("Chebyshev::parallel::CorrelationExpansionMoments will used %f GB\n", chebVecL.MemoryConsumptionInGB() + chebVecR.MemoryConsumptionInGB() );
      
  std::cout<<"Initialize chebVecL"<<std::endl;
  chebVecL=chebyshev::Vectors( chebMoms );
  std::cout<<"Initialize chebVecR"<<std::endl;
  chebVecR=chebyshev::Vectors( chebMoms );
  //}
  
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
	    if( j <10)
	      std::cout<<j<<" "<<Phi[j]<<std::endl;
	    
	    break;
	  case LOCAL_STATE:
 //	    i=(numStates-1); 
	    Phi[j] = ( (j==i) ? 1.0 : 0.0 ) ;
	    break;
	  default:
	    std::cerr<<" The state state is not identify, aborting running"<<std::endl;
	    std::exit(-1);
	  }
//	numStates= 1;
	//SELECT STATE TYPE
	std::cout<<"Computing with ID: "<<i+1<<" of "<<numStates<<" states" <<std::endl;

      //SELECT RUNNING TYPE
      /*if( use_sequential )
	chebyshev::sequential::CorrelationExpansionMoments(Phi, Phi, HAM, OPL, OPR, chebMoms);
	else*/
      chebyshev::parallel::TempDensityExpansionMoments(batchSize, Phi, HAM, OP, chebVecL, chebVecR, chebMoms);
    }
  
  //Fix the scaling of the moments
  const int NumMoms = chebMoms.HighestMomentNumber();
  const int NumTimes = chebMoms.HighestTimeNumber();
  for (int m = 0 ; m < NumMoms; m++)				  
    for (int n = m; n < NumTimes; n++)
      {
	double scal = 4.0/numStates;
	if( m == 0) scal*=0.5;
		
	const value_t tmp = scal * ( chebMoms(m, n) + std::conj(chebMoms(m, n)) ) / 2.0;
	chebMoms(m, n) = tmp;
	chebMoms(m, n) = std::conj(tmp);
      }
  
  return 0;
};
