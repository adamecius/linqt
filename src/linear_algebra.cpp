
#include "linear_algebra.hpp"

void linalg::axpy(const int dim, complex<double> a, const complex<double> *x, complex<double> *y)
{
	cblas_zaxpy(dim, &a, x, 1, y, 1);
}

void linalg::copy(const int dim, const complex<double> *x, complex<double> *y)
{
	cblas_zcopy(dim, x, 1, y, 1);
}


complex<double> linalg::vdot(const int dim, const complex<double> *x, complex<double> *y)
{
	complex<double> dotc;
	cblas_zdotc_sub(dim, x, 1, y, 1, &dotc);
	return dotc;
}


void linalg::batch_vdot(const int dim,const int batchSize,const complex<double>* leftb,const complex<double>* rightb,complex<double>* output)
{
	complex<double> alpha=1.0, beta=1.0; //This gives the quantity <L|R>* because there is no pure conjugation in MKL
	cblas_zgemm (CblasRowMajor, CblasNoTrans, CblasConjTrans,batchSize, batchSize , dim, &alpha, leftb,dim ,rightb,dim,&beta,output,batchSize);
	return ;
}


