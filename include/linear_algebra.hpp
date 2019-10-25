

#ifndef ALGEBRA_FUNCTION
#define ALGEBRA_FUNCTION
#include <complex>
using namespace std;


#define MKL_Complex16 complex<double>
#include "mkl.h"


namespace linalg
{

// VECTOR-VECTOR FUNCTIONS
void axpy(const int dim , complex<double> a,const complex<double>* x, complex<double>* y );

void copy(const int dim , const complex<double>* x, complex<double>* y );

complex<double> vdot(const int dim , const complex<double>* x, complex<double>* y );

}

#endif
