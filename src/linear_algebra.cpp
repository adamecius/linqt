
#include "linear_algebra.hpp"


void linalg::axpy(const int dim , complex<double> a,const complex<double>* x, complex<double>* y )
{
	cblas_zaxpy(dim, &a,x,1,y,1);
}

void linalg::copy(const int dim , const complex<double>* x, complex<double>* y )
{
	cblas_zcopy(dim,x,1,y,1);
}

complex<double> linalg::vdot(const int dim , const complex<double>* x, complex<double>* y )
{
	complex<double> dotc;
	cblas_zdotc_sub (dim ,x,1,y,1, &dotc);
	return dotc;
}
