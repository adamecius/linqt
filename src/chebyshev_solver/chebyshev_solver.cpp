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
    int batchSize = ( ( !std::string(getenv("BATCH_SIZE")).empty() ) ? atoi(getenv("BATCH_SIZE")):3);
    if( batchSize < 3){ std::cout<<"The minimum batch size is 3, so setting to 3"<<std::endl; batchSize=3;}
    if( batchSize > MAXSIZE){ std::cout<<"The batch larger than the number of moments is a waste of resources. Therefore set to :"<<MAXSIZE<<std::endl; batchSize=MAXSIZE;}
    std::cout<<"Using Bath Size of : "<<batchSize<<std::endl;
    std::cout<<"Allocating: "<<( (3*batchSize)*DIM + MAXSIZE*MAXSIZE)*sizeof(complex<double>)*1e-9<<"GB"<<std::endl;
    complex<double> *data = new complex<double>[(3*batchSize)* DIM + batchSize*batchSize];
    complex<double> **JL  = new complex<double> *[batchSize];
    complex<double> **JR  = new complex<double> *[batchSize];
    for (int b = 0; b < batchSize; b++)
    {
        JL[b] = &data[(b+0 * batchSize) * DIM];
        JR[b] = &data[(b+1 * batchSize) * DIM];
    }
    complex<double> *JV        = &data[2*batchSize*DIM];
    complex<double> *tmp_table = &data[3*batchSize*DIM];
    complex<double> *Jswap; //
    std::cout<<"memory allocated"<<std::endl;


    //INITIALIZE ITERATION
    vector<complex<double> > Phi(DIM);
    for (int i = 0; i < numStates; i++)
    {
        //while (stateFactory.CreateNewState())

        for (int i = 0; i < DIM; i++)
        {
            const complex<double> I(0, 1);
            Phi[i] = exp(I *2.0*M_PI* (double)rand() / (double)RAND_MAX)/sqrt(DIM);
        }

        //Start the chebyshev expansion of the correlations
        OPR.Multiply(1.0, &Phi[0], 0.0, JR[0]);
        HAM.Multiply(scalFactor, JR[0], 0.0, JR[1]);
        linalg::axpy(DIM, -shift,JR[0], JR[1]);
        for (int m1 = 0; m1 < cTable.Size_InDir(1); m1 += batchSize)
        {
	    for (int mR = 2; mR < batchSize; mR++)
            if (mR + m1 < cTable.Size_InDir(1))
            {
		linalg::copy(DIM, JR[mR-2], JR[mR]);
		HAM.Multiply(2.0 * scalFactor, JR[mR-1], -1.0, JR[mR]);
		linalg::axpy(DIM, -2.0*shift, JR[mR-1], JR[mR]);
	    }
            linalg::copy(DIM, &Phi[0], JL[0]);
            HAM.Multiply(scalFactor, JL[0], 0.0, JL[1]);
	    linalg::axpy(DIM, -shift, JL[0], JL[1]);
            for (int m0= 0; m0 < cTable.Size_InDir(0); m0 += batchSize)
            {
                for (int mL = 2; mL < batchSize; mL++)
                if (mL + m0 < cTable.Size_InDir(0))
                {
                   linalg::copy(DIM, JL[mL-2], JL[mL]);
                   HAM.Multiply(2.0 * scalFactor, JL[mL-1], -1.0, JL[mL]);
                   linalg::axpy(DIM, -2.0*shift,  JL[mL-1], JL[mL]);
                }
                //Compute the moment
		OPL.BatchMultiply(batchSize,1.0,*JL ,0.0, JV );
		linalg::batch_vdot(DIM,batchSize,*JR,JV,tmp_table ); //This actually gives <JR|JL>*

                for (int mR = 0; mR < batchSize; mR++)
                    for (int mL = 0; mL < batchSize; mL++)
                        if (mR+m1 < cTable.Size_InDir(1) && mL+m0 < cTable.Size_InDir(0))
			{
			   double scal=4;
			   if( mL+m0==0) scal*=0.5;
			   if( mR+m1==0) scal*=0.5;
			   cTable(mL+m0,mR+m1)+=( std::conj(tmp_table[mL*batchSize + mR])+tmp_table[mR*batchSize + mL] )*scal/2.0;
			}

                Jswap = JL[batchSize - 2];
                JL[batchSize - 2] = JL[0];
                JL[0] = Jswap;
                Jswap = JL[batchSize - 1];
                JL[batchSize - 1] = JL[1];
                JL[1] = Jswap;
            }
            Jswap = JR[batchSize - 2];
            JR[batchSize - 2] = JR[0];
            JR[0] = Jswap;
            Jswap = JR[batchSize - 1];
            JR[batchSize - 1] = JR[1];
            JR[1] = Jswap;
        }
    }
    delete[] JL;
    delete[] JR;
    delete[] data;
};
