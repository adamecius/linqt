
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

void chebyshev::CorrelationExpansionMoments(Factory &stateFactory, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR, chebyshev::MomTable &cTable)
{
    //Configure chebyshev
    const int MAXSIZE = cTable.maxSize();
    cTable.SetSystemSize(HAM.rank());
    const int
        DIM = HAM.rank(),
        MAX_MEM = cTable.Configure().maxMemory,
        VECSIZE = sizeof(complex<double>) * 2 * DIM; //The two is because one need two set of vectors

    int batchSize = (MAX_MEM + VECSIZE - 1) / VECSIZE;
    if (batchSize > MAXSIZE)
        batchSize = MAXSIZE;

    const double scalFactor = cTable.ScaleFactor();
    const double shift = cTable.EnergyShift();
    //Allocate the memory
    complex<double> *data = new complex<double>[2 * batchSize * DIM];
    complex<double> *JL = &data[0 * batchSize * DIM];
    complex<double> *JR = &data[1 * batchSize * DIM];
    vector<complex<double> > Phi(DIM);

    //INITIALIZE ITERATION
    while (stateFactory.CreateNewState())
    {
        Phi = stateFactory.getNewState();

        //Start the chebyshev expansion of the correlations
        OPR.Multiply(1.0, &Phi[0], 0.0, &JR[0]);
        HAM.Multiply(scalFactor, &JR[0], 0.0, &JR[1]);
        linalg::axpy(DIM, shift, &JR[0], &JR[1]);
        for (int m0 = 0; m0 < cTable.Size_InDir(0); m0 += batchSize)
        {
            for (int mR = m0 + 2; mR < m0 + batchSize; mR++)
                if (mR < cTable.Size_InDir(0))
                {
                    linalg::copy(DIM, &JR[mR - 2], &JR[mR]);
                    HAM.Multiply(2.0 * scalFactor, &JR[mR - 1], -1.0, &JR[mR]);
                    linalg::axpy(DIM, shift, &JR[mR - 1], &JR[mR]);
                }

            linalg::copy(DIM, &Phi[0], &JL[0]);
            HAM.Multiply(scalFactor, &JL[0], 0.0, &JL[1]);
            for (int m1 = 0; m1 < cTable.Size_InDir(1); m1 += batchSize)
            {
                for (int mL = m1 + 2; mL < m1 + batchSize; mL++)
                    if (mL < cTable.Size_InDir(1))
                    {
                        linalg::copy(DIM, &JL[mL - 2], &JL[mL]);
                        HAM.Multiply(scalFactor, &JL[mL - 1], -1.0, &JL[mL]);
                        linalg::axpy(DIM, shift, &JL[mL - 1], &JR[mL]);
                    }
                //Compute the moments
                for (int mR = 0; mR < m0 + batchSize; mR++)
                    for (int mL = 0; mL < m1 + batchSize; mL++)
                        if (mR < cTable.Size_InDir(0) && mL < cTable.Size_InDir(1))
                            cTable(mL, mR) += linalg::vdot(DIM, &JR[mL], &JL[mR]);
            }
        }
    }
    delete[] data;
};

