#include "linear_algebra.hpp"
#include <eigen-3.4.0/Eigen/Core>

void linalg::scal(const complex<double>& a, vector< complex<double> >& x)
{
	cblas_zscal(x.size(), &a, &x[0], 1);
	return ;
}


void linalg::axpy(const complex<double>& a, const vector< complex<double> >& x, vector< complex<double> >& y)
{
	assert( x.size() == y.size() );
	cblas_zaxpy(x.size(), &a, &x[0], 1, &y[0], 1);
	return ;
}

void linalg::copy(const vector< complex<double> >&x,vector< complex<double> >& y)
{
	assert( x.size() == y.size() );
	cblas_zcopy(x.size(), &x[0], 1, &y[0], 1);
}

complex<double> linalg::vdot(const vector< complex<double> >& x,const vector< complex<double> >& y)
{
	assert( x.size() == y.size() );
	complex<double> dotc;
	cblas_zdotc_sub(x.size(), &x[0], 1, &y[0], 1, &dotc);
	return dotc;
}

double linalg::nrm2(const vector< complex<double> >& x)
{
	return  cblas_dznrm2 (x.size(),&x[0],1);
};


void linalg::scal(const int dim, complex<double> a, complex<double> *x)
{
	cblas_zscal(dim, &a, x, 1);
	return ;
}

void linalg::axpy(const int dim, complex<double> a, const complex<double> *x, complex<double> *y)
{
	cblas_zaxpy(dim, &a, x, 1, y, 1);
	return ;
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


double linalg::nrm2(const int dim, const complex<double> *x)
{
	return  cblas_dznrm2 (dim,x,1);
};

void linalg::batch_vdot(const int dim,const int batchSize,const complex<double>* leftb,const complex<double>* rightb,complex<double>* output)
{
	complex<double> alpha=1.0, beta=1.0; //This gives the quantity <L|R>* because there is no pure conjugation in MKL
	cblas_zgemm (CblasRowMajor, CblasNoTrans, CblasConjTrans,batchSize, batchSize , dim, &alpha, leftb,dim ,rightb,dim,&beta,output,batchSize);

	const long tot_dim = batchSize*batchSize;
	#pragma omp parallel for 
	for(long int j =0 ; j < tot_dim ; j++)
		output[j]= std::conj(output[j]);
	return ;
}



void linalg::extract_segment(const vector< complex<double> >&x, size_t size_x, size_t start_x,  vector< complex<double> >& y, size_t size_y ){//size_x >> size_y

#pragma omp parallel for 
  for(size_t i = 0;  i < size_y ; i++)//&& (start_x + i < size_x)  !!!
    y[i] = x[start_x + i];

};

void linalg::introduce_segment(const vector< complex<double> >&x, size_t size_x, vector< complex<double> >& y, size_t size_y, size_t start_y ){//size_y >> size_y


#pragma omp parallel for 
  for(size_t i = 0; i< size_x; i++)
    y[start_y + i] = x[i];

};  


void linalg::orthogonalize(SparseMatrixType& S, vector< complex<double> >& orthogonalized, vector< complex<double> >& original){


  int DIM = S.numRows(), NNZ=S.vals()->size();


  int *rowsPtr = S.rows()->data();
  int *colsPtr = S.cols()->data();
  std::complex<double> *values = S.vals()->data();


  Eigen::Map<Eigen::Vector<std::complex<double>, -1>>
    eig_original(original.data(), DIM),
    eig_orthogonalized(orthogonalized.data(), DIM);



  Eigen::Map<Eigen::SparseMatrix<complex<double>, Eigen::RowMajor> > eigen_S( DIM, DIM, NNZ, rowsPtr, colsPtr, values);

  
  Eigen::BiCGSTAB<Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> > solver;
  solver.setTolerance(0.0001); 
  solver.setMaxIterations(2000); 
  solver.compute(eigen_S);
  eig_orthogonalized = solver.solve(eig_original);


  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "  max#iterations:" << solver.maxIterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;
  std::cout << "  tolerance :    " << solver.tolerance()      << std::endl;
  std::cout<<  "Vector norm :    " <<eig_original.norm()<<std::endl<<std::endl;    
};
