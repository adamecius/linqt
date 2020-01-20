// Used for OPENMP functions
#include "chebyshev_solver.hpp"
#include <omp.h>

bool chebyshev::GetBatchSize(int& batchSize)
{
	batchSize = 3;
	if(!getenv("BATCH_SIZE") || atoi(getenv("BATCH_SIZE")) <= 0) //evaluation left-to-right
		std::cout<<"\nEnviroment variable BATCH_SIZE not set or invalid.\nSet to your custom N value throught the command export BATCH_SIZE N. Using sequential"<<std::endl;
	else
	{
		batchSize = atoi(getenv("BATCH_SIZE"));
		return true;
	}	
	return false;
}

void chebyshev::MomTable::saveIn(std::string filename)
{
    const double bandWidth = 2/(double)conf_->scaleFactor;
    const double bandCenter = conf_->shift/(double)conf_->scaleFactor;

    typedef std::numeric_limits<double> dbl;

    ofstream outputfile(filename.c_str());
    outputfile.precision(dbl::digits10);
    outputfile << systSize_ << " " << bandWidth << "  " << bandCenter << std::endl;
    //Print the number of moments for all directions in a line
    for (vector<int>::iterator it = size_.begin();
         it != size_.end(); it++)
        outputfile << *it << " ";
    outputfile << std::endl;

    for (vector<complex<double> >::iterator it = data_.begin();
         it != data_.end(); it++)
        outputfile << (*it).real() << " " << (*it).imag() << std::endl;
    outputfile.close();
};

chebyshev::Moments2D::Moments2D( std::string momfilename )
{
	//Check if the input_momfile have the right extension 
	std::size_t ext_pos = string( momfilename ).find(".chebmom2D"); 
	if( ext_pos == string::npos )
	{ std::cerr<<"The first argument does not seem to be a valid .chebmom2D file"<<std::endl; assert(false);}

	//if it does, use it to get the extension
	system_label = momfilename.substr(0,ext_pos);


	//and then try to open the file
	std::ifstream momfile( momfilename.c_str() );
	assert( momfile.is_open() );

	//if succesful, read the header
	momfile>>this->system_size>>this->band_width>>this->band_center; //in the file what you have is the bandwidth
	momfile>>this->numMoms[0]>>this->numMoms[1];
		
	//create the moment array and read the data
	mu = vector<complex<double> >(numMoms[1]*numMoms[0], 0.0);	
	double rmu,imu;
	for( int m0 = 0 ; m0 < numMoms[0] ; m0++)
	for( int m1 = 0 ; m1 < numMoms[1] ; m1++)
	{ 
		momfile>>rmu>>imu;
		this->operator()(m0,m1) = std::complex<double>(rmu,imu);
	}
	momfile.close();
};


void chebyshev::DensityMoments( std::vector< std::complex<double> >& PhiL,
								std::vector< std::complex<double> >& PhiR,
								SparseMatrixType &HAM, 
								SparseMatrixType &OPL, 
								const int numMoms,
								const double scalFactor,
								const double shift,
								std::complex<double>* mu)
{ 
    const int DIM = HAM.rank();
	std::vector< std::complex<double> > J0(DIM), J1(DIM), JOPL(DIM);
	
	linalg::copy(DIM, &PhiL[0], &J0[0]);
	HAM.Multiply(scalFactor, &J0[0], 0.0, &J1[0]);
	linalg::axpy(DIM, -shift,&J0[0],&J1[0]);
	for(int m=0; m < numMoms ; m++)
	{
		OPL.Multiply(1.0,&J0[0] ,0.0, &JOPL[0] );
		mu[m] += linalg::vdot(DIM,&JOPL[0], &PhiR[0] ); //This actually gives <JR|JL>*
		HAM.Multiply(2.0*scalFactor, &J1[0], -1.0,&J0[0]);
		linalg::axpy(DIM,-2.0*shift, &J1[0], &J0[0]);
		J0.swap(J1);
	}
};


void chebyshev::Vectors( 					 
						  std::vector< std::complex<double> >& J0,
						  std::vector< std::complex<double> >& J1,
						 SparseMatrixType &HAM,
						 const int numMoms,
						 const double scalFactor,
						 const double shift,
						 std::vector< std::vector< std::complex<double> > >& JCheb)
{
	assert( JCheb.size()  >= numMoms );
	const int DIM = HAM.rank();
	std::complex<double> *JCheb0, *JCheb1, *JChebm ;

	for(int m = 0 ; m < numMoms; m++ )
	{
		JChebm=  &((JCheb[m])[0]);	
		switch ( m )
		{
			case 0:
				JCheb0=&(J0[0]);	
				JCheb1=&(J1[0]);	
				break;
			case 1:
				JCheb0=&(J1[0]);	
				JCheb1=&((JCheb[m-1])[0]);	
				break;
			default:
				JCheb0=&((JCheb[m-2])[0]);	
				JCheb1=&((JCheb[m-1])[0]);	
		}
		linalg::copy(DIM, JCheb0, JChebm);
		HAM.Multiply(2.0*scalFactor, JCheb1, -1.0,JChebm);
		linalg::axpy(DIM,-2.0*shift, JCheb1, JChebm);
	}
};

void chebyshev::sequential::CorrelationExpansionMoments(const int numMoms0,const int numMoms1,
											const std::vector< std::complex<double> >& Phi,
											SparseMatrixType &HAM, 
											SparseMatrixType &OPL, 
											SparseMatrixType &OPR, 
											chebyshev::MomTable &cTable)
{
    const int DIM = HAM.rank();
    const double scalFactor = cTable.ScaleFactor();
    const double shift = cTable.EnergyShift();
	std::vector< std::complex<double> > JR0(DIM),JR1(DIM),JL0(DIM),JL1(DIM),JOL(DIM);

    //Start the chebyshev expansion of the correlations
	OPR.Multiply(1.0, &Phi[0], 0.0, &JR0[0]);
    HAM.Multiply(scalFactor,&JR0[0], 0.0, &JR1[0]);
    linalg::axpy(DIM,-shift,&JR0[0], &JR1[0]);
    for (int m1 = 0; m1 < numMoms1; m1++ )
    {
		linalg::copy(DIM, &Phi[0], &JL0[0]);
		HAM.Multiply(scalFactor, &JL0[0], 0.0, &JL1[0]);
		linalg::axpy(DIM, -shift,&JL0[0],&JL1[0]);
		for (int m0= 0; m0 < numMoms0; m0++)
		{
			OPL.Multiply(1.0,&JL0[0] ,0.0, &JOL[0] );
			cTable(m0,m1) += linalg::vdot(DIM,&JOL[0], &JR0[0] ); //This actually gives <JR|JL>*
			
			HAM.Multiply(2.0 * scalFactor, &JL1[0], -1.0, &JL0[0]);
			linalg::axpy(DIM,-2.0*shift, &JL1[0], &JL0[0]);
			JL0.swap(JL1);
		}
		HAM.Multiply(2.0 * scalFactor, &JR1[0], -1.0, &JR0[0]);
		linalg::axpy(DIM,-2.0*shift, &JR1[0], &JR0[0]);
		JR0.swap(JR1);
	}
};


void chebyshev::CorrelationExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR, chebyshev::MomTable &cTable)
{
	//Allocate the memory
    string NUM_THREADS ="Default"; 	if(getenv("OMP_NUM_THREADS")) NUM_THREADS = getenv("OMP_NUM_THREADS");
	int batchSize=3;
    bool use_sequential = !( chebyshev::GetBatchSize(batchSize) );

	const double scalFactor = cTable.ScaleFactor();
    const double shift = cTable.EnergyShift();
    const int DIM = HAM.rank();

    cTable.SetSystemSize(HAM.rank());
    const int NumMomsL = cTable.Size_InDir(0);
    const int NumMomsR = cTable.Size_InDir(1);

    const double total_memory = ( ((double)batchSize + 5 )*(double)DIM +NumMomsL*(double)(NumMomsR) )*(double)sizeof(complex<double>)/pow(2.0,30.0);
    std::cout<<"Allocating: "<<total_memory<<"GB"<<std::endl;

	if( use_sequential )
	{
		std::cout<<"USING SEQUENTIAL IMPLEMENTATION"<<std::endl;
		std::vector< std::complex<double> > Phi(DIM); 	//States Vectors
		for (int i = 0; i < numStates; i++)
		{
			//construct a normalized state for the left side
			for (int j = 0; j < DIM; j++)
			{
				const complex<double> I(0, 1);
				Phi[j] = exp(I*2.0*M_PI* (double)rand() / (double)RAND_MAX)/sqrt(DIM);
			}
			
			chebyshev::sequential::CorrelationExpansionMoments(	NumMomsL,NumMomsR, Phi, HAM, OPL, OPR, cTable);
		}
	}
	else
	{
		std::cout<<"LARGE MEMORY CALCULATION "<<std::endl
				 <<"Using batch of "<<batchSize<<" in  "<<NUM_THREADS<<" threads"<<std::endl;

		std::vector< std::vector< std::complex<double> > > JR(batchSize);
		for(int i = 0 ; i < batchSize; i++ )
			JR[i] = std::vector< std::complex<double> >( DIM ); 
		int end = batchSize-1;
		
		//First we go throguh all the states
		std::vector< std::complex<double> > PhiR(DIM),PhiL(DIM); 	//States Vectors
		for (int i = 0; i < numStates; i++)
		{
			std::cout<<"Computing state: "<<i<<" out of "<<numStates<<std::endl;
			//construct a normalized state for the left side
			for (int j = 0; j < DIM; j++)
			{
				const complex<double> I(0, 1);
				PhiL[j] = exp(I*2.0*M_PI* (double)rand() / (double)RAND_MAX)/sqrt(DIM);
			}

			//construct the right side state
			OPR.Multiply(1.0, &PhiL[0], 0.0, &PhiR[0]);

			//COMPUTE BATCH OF LEFT VECTORS
			for(int  mR = 0 ; mR <  NumMomsR ; mR+=batchSize)
			{
				switch ( mR )
				{
					case 0:
						linalg::copy(DIM, &PhiR[0], &(JR[mR][0]) );
						chebyshev::DensityMoments( PhiL , JR[mR], HAM, OPL,NumMomsL,scalFactor, shift, &cTable(0, mR  ) );
						mR++;

					case 1:
						HAM.Multiply(scalFactor, &(JR[mR-1][0]),0.0,&(JR[mR][0]) );
						linalg::axpy(DIM, -shift,&(JR[mR-1][0]),    &(JR[mR][0]));
						chebyshev::DensityMoments( PhiL , JR[mR], HAM, OPL,NumMomsL,scalFactor, shift, &cTable(0, mR  ) );
						mR++;
						
						//At this point we move the data to the end of the array because is going to be necessary
						//for iterations
						linalg::copy(DIM, &(JR[mR-2][0]),&(JR[end-1][0]) );//We move the initial vectors to the end position - 1.
						linalg::copy(DIM, &(JR[mR-1][0]),&(JR[end  ][0]) );//We move the initial vectors to the end position.

					default:
						chebyshev::Vectors( JR[end-1],JR[end] ,HAM, batchSize, scalFactor, shift,JR);
					break;
				}
				for(int m0= 0 ; m0 < batchSize ; m0++)
				if( mR + m0 < NumMomsR )
					chebyshev::DensityMoments( PhiL , JR[m0], HAM, OPL,NumMomsL,scalFactor, shift, &cTable(0, mR + m0 ) );
			}
		}	
	}

	for (int mL = 0 ; mL < NumMomsL; mL++)				  
	for (int mR = mL; mR < NumMomsR; mR++)
	{
		double scal=4.0/numStates;
		if( mL==0) scal*=0.5;
		if( mR==0) scal*=0.5;

		const std::complex<double> tmp = scal*( cTable(mL,mR) + std::conj(cTable(mR,mL)) )	/2.0;
		cTable(mL,mR)= tmp;
		cTable(mR,mL)= std::conj(tmp);
		
	}

};

void chebyshev::LocalCorrelationExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR, chebyshev::MomTable &cTable)
{
	//Allocate the memory
    string NUM_THREADS ="Default"; 	if(getenv("OMP_NUM_THREADS")) NUM_THREADS = getenv("OMP_NUM_THREADS");
	int batchSize=3;
    bool use_sequential = !( chebyshev::GetBatchSize(batchSize) );

	const double scalFactor = cTable.ScaleFactor();
    const double shift = cTable.EnergyShift();
    const int DIM = HAM.rank();

    cTable.SetSystemSize(HAM.rank());
    const int NumMomsL = cTable.Size_InDir(0);
    const int NumMomsR = cTable.Size_InDir(1);

    const double total_memory = ( ((double)batchSize + 5 )*(double)DIM +NumMomsL*(double)(NumMomsR) )*(double)sizeof(complex<double>)/pow(2.0,30.0);
    std::cout<<"Allocating: "<<total_memory<<"GB"<<std::endl;

	if( use_sequential )
	{
		std::cout<<"USING SEQUENTIAL IMPLEMENTATION"<<std::endl;
		std::vector< std::complex<double> > Phi(DIM); 	//States Vectors
		for (int i = 0; i < numStates; i++)
		{
			//construct a normalized state for the left side
			for (int j = 0; j < DIM; j++)
				Phi[j] = ( (j==i) ? 1.0 : 0.0 ) ;
			
			chebyshev::sequential::CorrelationExpansionMoments(	NumMomsL,NumMomsR, Phi, HAM, OPL, OPR, cTable);
		}
	}
	else
	{
		std::cout<<"LARGE MEMORY CALCULATION "<<std::endl
				 <<"Using batch of "<<batchSize<<" in  "<<NUM_THREADS<<" threads"<<std::endl;

		std::vector< std::vector< std::complex<double> > > JR(batchSize);
		for(int i = 0 ; i < batchSize; i++ )
			JR[i] = std::vector< std::complex<double> >( DIM ); 
		int end = batchSize-1;
		
		//First we go throguh all the states
		std::vector< std::complex<double> > PhiR(DIM),PhiL(DIM); 	//States Vectors
		for (int i = 0; i < numStates; i++)
		{
			std::cout<<"Computing state: "<<i<<" out of "<<numStates<<std::endl;
			for (int j = 0; j < DIM; j++)
				PhiL[j] = ( (j==i) ? 1.0 : 0.0 ) ;

			//construct the right side state
			OPR.Multiply(1.0, &PhiL[0], 0.0, &PhiR[0]);

			//COMPUTE BATCH OF LEFT VECTORS
			for(int  mR = 0 ; mR <  NumMomsR ; mR+=batchSize)
			{
				switch ( mR )
				{
					case 0:
						linalg::copy(DIM, &PhiR[0], &(JR[mR][0]) );
						chebyshev::DensityMoments( PhiL , JR[mR], HAM, OPL,NumMomsL,scalFactor, shift, &cTable(0, mR  ) );
						mR++;

					case 1:
						HAM.Multiply(scalFactor, &(JR[mR-1][0]),0.0,&(JR[mR][0]) );
						linalg::axpy(DIM, -shift,&(JR[mR-1][0]),    &(JR[mR][0]));
						chebyshev::DensityMoments( PhiL , JR[mR], HAM, OPL,NumMomsL,scalFactor, shift, &cTable(0, mR  ) );
						mR++;
						
						//At this point we move the data to the end of the array because is going to be necessary
						//for iterations
						linalg::copy(DIM, &(JR[mR-2][0]),&(JR[end-1][0]) );//We move the initial vectors to the end position - 1.
						linalg::copy(DIM, &(JR[mR-1][0]),&(JR[end  ][0]) );//We move the initial vectors to the end position.

					default:
						chebyshev::Vectors( JR[end-1],JR[end] ,HAM, batchSize, scalFactor, shift,JR);
					break;
				}
				for(int m0= 0 ; m0 < batchSize ; m0++)
				if( mR + m0 < NumMomsR )
					chebyshev::DensityMoments( PhiL , JR[m0], HAM, OPL,NumMomsL,scalFactor, shift, &cTable(0, mR + m0 ) );
			}
		}	
	}

	for (int mL = 0 ; mL < NumMomsL; mL++)				  
	for (int mR = mL; mR < NumMomsR; mR++)
	{
		double scal=4.0/numStates;
		if( mL==0) scal*=0.5;
		if( mR==0) scal*=0.5;

		const std::complex<double> tmp = scal*( cTable(mL,mR) + std::conj(cTable(mR,mL)) )	/2.0;
		cTable(mL,mR)= tmp;
		cTable(mR,mL)= std::conj(tmp);
		
	}
	
}








/*
void chebyshev::CorrelationExpansionMoments(int numStates, 
											SparseMatrixType &HAM, 
											SparseMatrixType &OPL, 
											SparseMatrixType &OPR, 
											chebyshev::MomTable &cTable)
{
    //Configure chebyshev
    const int MAXSIZE = cTable.maxSize();
    cTable.SetSystemSize(HAM.rank());
    const int DIM = HAM.rank();

    const double scalFactor = cTable.ScaleFactor();
    const double shift = cTable.EnergyShift();

    //Allocate the memory
	if(!getenv("BATCH_SIZE"))
		std::cout<<"\nEnviroment variable BATCH_SIZE not set, using BATCH_SIZE=1.\nSet to your custom N value throught the command export BATCH_SIZE N"<<std::endl;
	else if ( atoi(getenv("BATCH_SIZE"))< 0 )
		std::cout<<"\Enviroment variable BATCH_SIZE is negative, set to 1, the minimum allowed. \nSet to your custom N value throught the command export BATCH_SIZE N"<<std::endl;
		
		
    int batchSize = ( ( getenv("BATCH_SIZE") ) ? atoi(getenv("BATCH_SIZE")):1);
    if( batchSize < 0){ std::cout<<"The minimum batch size is 1, so setting to 1"<<std::endl; batchSize=1;}
    if( batchSize > MAXSIZE){ std::cout<<"The batch larger than the number of moments is a waste of resources. Therefore set to :"<<MAXSIZE<<std::endl; batchSize=MAXSIZE;}
    std::cout<<"Using Bath Size of : "<<batchSize<<std::endl;


    const double total_memory = ( (3*(double)batchSize + 5 )*(double)DIM + (double)batchSize*(double)batchSize)*(double)sizeof(complex<double>)/pow(2.0,30.0);
    std::cout<<"Allocating: "<<total_memory<<"GB"<<std::endl;
    std::vector< std::complex<double> > JR0(DIM), JR1(DIM), JL0(DIM), JL1(DIM), Phi(DIM);
    std::vector< std::complex<double> > tmp_table(batchSize*batchSize);	
    std::cout<<"MEMORY ALLOCATED "<<std::endl;

    //Set the correct address for the pointer to pointers
    vector< complex<double> > data( (long int)3*(long int) batchSize * (long int) DIM );
    complex<double>**JL  = new complex<double>*[batchSize];
    complex<double>**JR  = new complex<double>*[batchSize];
    complex<double>**JV  = new complex<double>*[batchSize];
    for (int b = 0; b < batchSize; b++)
    {
       JL[b] = &data[(b + 0*batchSize )*DIM];
       JR[b] = &data[(b + 1*batchSize )*DIM];
       JV[b] = &data[(b + 2*batchSize )*DIM];//dimension batchSize*DIM
    }

    //INITIALIZE ITERATION
    for (int i = 0; i < numStates; i++)
    {
        //while (stateFactory.CreateNewState())

        for (int i = 0; i < DIM; i++)
        {
            const complex<double> I(0, 1);
            Phi[i] = exp(I*2.0*M_PI* (double)rand() / (double)RAND_MAX)/sqrt(DIM);
        }

        //Start the chebyshev expansion of the correlations
        OPR.Multiply(1.0, &Phi[0], 0.0, &JR0[0]);
        HAM.Multiply(scalFactor,&JR0[0], 0.0, &JR1[0]);
        linalg::axpy(DIM,-shift,&JR0[0], &JR1[0]);
        for (int m1 = 0; m1 < cTable.Size_InDir(1); m1 += batchSize)
        {
			for (int mR = 0; mR < batchSize; mR++)
			if (mR + m1 < cTable.Size_InDir(1))
			{
				linalg::copy(DIM, &JR0[0], JR[mR]);
				HAM.Multiply(2.0 * scalFactor, &JR1[0], -1.0, &JR0[0]);
				linalg::axpy(DIM,-2.0*shift, &JR1[0], &JR0[0]);
				JR0.swap(JR1);
			}

            //Start the inner ctable matrix loop
            linalg::copy(DIM, &Phi[0], &JL0[0]);
            HAM.Multiply(scalFactor, &JL0[0], 0.0, &JL1[0]);
            linalg::axpy(DIM, -shift,&JL0[0],&JL1[0]);
            for (int m0= 0; m0 < cTable.Size_InDir(0); m0 += batchSize)
            {
				for (int mL = 0; mL < batchSize; mL++)
				if (mL + m0 < cTable.Size_InDir(0))
				{
					linalg::copy(DIM, &JL0[0], JL[mL]);
					HAM.Multiply(2.0 * scalFactor, &JL1[0], -1.0, &JL0[0]);
					linalg::axpy(DIM,-2.0*shift, &JL1[0], &JL0[0]);
					JL0.swap(JL1);
				}
				//Compute the moment (HARD PART)
//				for (int mL = 0; mL < batchSize; mL++)
//					OPL.Multiply(1.0,JL[mL] ,0.0, JV[mL] );				
//				for (int mR = 0; mR < batchSize; mR++)
//				for (int mL = 0; mL < batchSize; mL++)
//					tmp_table[mR*batchSize + mL]=linalg::vdot(DIM,JV[mL],JR[mR]); //This actually gives <JR|JL>*
				std::cout<<"Computing batch ("<<m0<<","<<m1<<")"<<std::endl;

				OPL.BatchMultiply(batchSize,1.0,*JL ,0.0, *JV );
				linalg::batch_vdot(DIM,batchSize,*JR,*JV,&tmp_table[0] ); //This actually gives <JR|JL>*			
				for (int mR = 0; mR < batchSize; mR++)				  // therefore, the tmp_table appears as conjugaated for the right order
				for (int mL = 0; mL < batchSize; mL++)
				if (mR+m1 < cTable.Size_InDir(1) && mL+m0 < cTable.Size_InDir(0))
				{
					double scal=4;
					if( mL+m0==0) scal*=0.5;
					if( mR+m1==0) scal*=0.5;
					cTable(mL+m0,mR+m1)+= scal*( std::conj(tmp_table[mR*batchSize + mL] ) + tmp_table[mL*batchSize + mR]  )/2.0;
				}
            }
        }
    }
    std::cout<<"Release batch vector's memory"<<std::endl;
    delete[] JL; JL=0;
    delete[] JR; JR=0;
    delete[] JV; JV=0;

};
*/

