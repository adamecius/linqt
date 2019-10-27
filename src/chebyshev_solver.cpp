
// Used for OPENMP functions
#include "chebyshev_solver.hpp"

void chebyshev::MomTable::saveIn(std::string filename)
{
    const double bandWidth = 2 * conf_->scaleFactor;
    const double bandCenter = -conf_->shift * conf_->scaleFactor;

    typedef std::numeric_limits<double> dbl;

    ofstream outputfile(filename.c_str());
    outputfile.precision(dbl::digits10);
    outputfile << systSize_ << " " << bandWidth << "  " << bandCenter << std::endl;
    //Print the number of moments for all directions in a line
    for (vector<int>::iterator it = size_.begin();
         it != size_.end(); it++)
        outputfile << *it << " " << std::endl;
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
    const int
        DIM = HAM.rank(),
        toGB = 1048576,
        VECSIZE = sizeof(complex<double>) * 2 * DIM / toGB; //The two is because one need two set of vectors
    int MAX_MEM = 40;
    int batchSize = (MAX_MEM + VECSIZE - 1) / VECSIZE;
    if (batchSize > MAXSIZE)
        batchSize = MAXSIZE;
    std::cout << "The dimension of my complex vector is: " << DIM << std::endl;
    std::cout << "Which will use  dimension of my complex vector is: " << VECSIZE << "GB " << std::endl;
    std::cout << "Because my Memory Size is =" << MAX_MEM << "I Can use a batch of " << batchSize << " vectors." << std::endl;

    const double scalFactor = cTable.ScaleFactor();
    const double shift = cTable.EnergyShift();
    //Allocate the memory
    complex<double> *data = new complex<double>[2 * batchSize * DIM];
    complex<double> **JL = new complex<double> *[batchSize];
    complex<double> **JR = new complex<double> *[batchSize];
    complex<double> *Jswap; //
    for (int b = 0; b < batchSize; b++)
    {
        JL[b] = &data[(b + 0 * batchSize) * DIM];
        JR[b] = &data[(b + 1 * batchSize) * DIM];
    }
    vector<complex<double> > Phi(DIM);

    //INITIALIZE ITERATION
    for (int i = 0; i < numStates; i++)
    {
        //while (stateFactory.CreateNewState())

        for (int i = 0; i < DIM; i++)
        {
            const complex<double> I(0, 1);
            Phi[i] = exp(I * (double)rand() / (double)RAND_MAX);
        }

        //Start the chebyshev expansion of the correlations
        OPR.Multiply(1.0, &Phi[0], 0.0, JR[0]);
        HAM.Multiply(scalFactor, JR[0], 0.0, JR[1]);
        linalg::axpy(DIM, shift, JR[0], JR[1]);
        for (int m0 = 0; m0 < cTable.Size_InDir(0); m0 += batchSize)
        {
            for (int mR = 2; mR < batchSize; mR++)
                if (mR + m0 < cTable.Size_InDir(0))
                {
                    linalg::copy(DIM, JR[mR - 2], JR[mR]);
                    HAM.Multiply(2.0 * scalFactor, JR[mR - 1], -1.0, JR[mR]);
                    linalg::axpy(DIM, shift, JR[mR - 1], JR[mR]);
                }

            linalg::copy(DIM, &Phi[0], JL[0]);
            HAM.Multiply(scalFactor, JL[0], 0.0, JL[1]);
            for (int m1 = 0; m1 < cTable.Size_InDir(1); m1 += batchSize)
            {
                for (int mL = 2; mL < batchSize; mL++)
                    if (mL + m1 < cTable.Size_InDir(1))
                    {
                        linalg::copy(DIM, JL[mL - 2], JL[mL]);
                        HAM.Multiply(2.0*scalFactor, JL[mL - 1], -1.0, JL[mL]);
                        linalg::axpy(DIM, shift, JL[mL - 1], JR[mL]);
                    }
                //Compute the moments
                for (int mR = 0; mR < batchSize; mR++)
                    for (int mL = 0; mL < batchSize; mL++)
                        if (mR < cTable.Size_InDir(0) && mL < cTable.Size_InDir(1))
                        {
                  //          std::cout<<"Computing moments: "<<mL+m1<<","<<mR+m0<<std::endl;
                            cTable(mL + m1, mR + m0) += linalg::vdot(DIM, JR[mL], JL[mR]);
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
