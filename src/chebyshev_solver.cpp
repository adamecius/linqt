// Used for OPENMP functions
#include "chebyshev_solver.hpp"

void chebyshev::MomTable::saveIn(std::string filename)
{
    const double bandWidth = 2/(double)conf_->scaleFactor;
    const double bandCenter = -conf_->shift/(double)conf_->scaleFactor;

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

void chebyshev::CorrelationExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR, chebyshev::MomTable &cTable)
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
