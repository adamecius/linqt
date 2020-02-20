#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX

#include <Eigen/Sparse>
#include <complex>

typedef Eigen::SparseMatrix< std::complex<double>, Eigen::RowMajor > SparseMatrix_t; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet< std::complex<double> > Triplet_t;


#endif