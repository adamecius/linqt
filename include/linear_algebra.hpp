

#ifndef ALGEBRA_FUNCTION
#define ALGEBRA_FUNCTION
#include <complex>
#include <vector>
#include <cassert>
#include "sparse_matrix.hpp"




namespace linalg
{

// VECTOR-STD::VECTOR FUNCTIONS
void scal(const int dim, std::complex<double> a, std::complex<double> *x);

void scal(const std::complex<double>& a, std::vector< std::complex<double> >& x);

void axpy(const int dim, std::complex<double> a, const std::complex<double> *x, std::complex<double> *y);

void axpy(const std::complex<double>& a, const std::vector< std::complex<double> >& x, std::vector< std::complex<double> >& y);

void copy(const int dim, const std::complex<double> *x, std::complex<double> *y);

void copy(const std::vector< std::complex<double> >&x,std::vector< std::complex<double> >& y);

std::complex<double> vdot(const std::vector< std::complex<double> >& x,const std::vector< std::complex<double> >& y);

std::complex<double> vdot(const int dim, const std::complex<double> *x, std::complex<double> *y);

double nrm2(const std::vector< std::complex<double> >& x);

double nrm2(const int dim, const std::complex<double> *x);

void batch_vdot(const int dim, const int batchSize, const std::complex<double>* leftb, const std::complex<double>* rightb, std::complex<double>* output);

void extract_segment( const std::vector< std::complex<double> >&x, size_t start_x, std::vector< std::complex<double> >& y);//size_x >> size_y
   
void introduce_segment( const std::vector< std::complex<double> >&x, std::vector< std::complex<double> >& y, size_t start_y );//size_y >> size_y

void orthogonalize(SparseMatrixType &,  const std::vector< std::complex<double> >&,  std::vector< std::complex<double> >& );
  
} // namespace linalg

#endif
