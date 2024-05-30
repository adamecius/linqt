#include "linear_algebra.hpp"
#include <eigen3/Eigen/Core>
#include<eigen3/Eigen/Sparse>
#include<eigen3/Eigen/IterativeLinearSolvers>

void linalg::scal(const std::complex<double>& a, std::vector< std::complex<double> >& x)
{
        Eigen::Map<Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x.data(), x.size());
            
       	eig_x*=a;
	
	return ;
}


void linalg::axpy(const std::complex<double>& a, const std::vector< std::complex<double> >& x, std::vector< std::complex<double> >& y)
{
	assert( x.size() == y.size() );
	
        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>> eig_x(x.data(), x.size());
	Eigen::Map<Eigen::Vector<std::complex<double>, -1>> eig_y(y.data(), y.size());
        
	eig_y = a * eig_x + eig_y;
	
	return ;
}

void linalg::copy(const std::vector< std::complex<double> >&x, std::vector< std::complex<double> >& y)
{
	assert( x.size() == y.size() );

        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>> eig_x(x.data(), x.size());
	Eigen::Map<Eigen::Vector<std::complex<double>, -1>> eig_y(y.data(), y.size());

	eig_y = eig_x;
	
}

std::complex<double> linalg::vdot(const std::vector< std::complex<double> >& x,const std::vector< std::complex<double> >& y)
{
	assert( x.size() == y.size() );
	
        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x.data(), x.size()),
	  eig_y(y.data(), y.size());

	return eig_x.dot(eig_y);
}

double linalg::nrm2(const std::vector< std::complex<double> >& x)
{
        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x.data(), x.size());

	return eig_x.squaredNorm();
};


void linalg::scal(const int dim, std::complex<double> a, std::complex<double> *x)
{
        Eigen::Map<Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x, dim);

       eig_x*=a;

	return ;
}

void linalg::axpy(const int dim, std::complex<double> a, const std::complex<double> *x, std::complex<double> *y)
{
        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x, dim);
	Eigen::Map<Eigen::Vector<std::complex<double>, -1>>
	  eig_y(y, dim);

	eig_y=a*eig_x+eig_y;
	
	return ;
}

void linalg::copy(const int dim, const std::complex<double> *x, std::complex<double> *y)
{
        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x, dim);

	Eigen::Map<Eigen::Vector<std::complex<double>, -1>>
	  eig_y(y, dim);

	eig_y = eig_x;
	
}

std::complex<double> linalg::vdot(const int dim, const std::complex<double> *x, std::complex<double> *y)
{	
        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x, dim);
	Eigen::Map<Eigen::Vector<std::complex<double>, -1>>
	  eig_y(y, dim);

	return eig_x.dot(eig_y);
}


double linalg::nrm2(const int dim, const std::complex<double> *x)
{
        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x, dim);

	return eig_x.squaredNorm();

};

void linalg::batch_vdot(const int dim, const int batchSize, const std::complex<double>* leftb, const std::complex<double>* rightb, std::complex<double>* output)
{
        Eigen::Map<const Eigen::MatrixXcd>
	   eig_leftb(leftb, dim, batchSize),
	   eig_rightb(rightb, dim, batchSize);
 
        Eigen::Map<Eigen::MatrixXcd> eig_output(output, batchSize, batchSize);

       	eig_output = eig_leftb.adjoint() * eig_rightb;

	return ;
}




void linalg::extract_segment(const std::vector< std::complex<double> >&x, size_t start_x,  std::vector< std::complex<double> >& y){//size_x >> size_y

#pragma omp parallel for 
  for(size_t i = 0;  i < y.size() ; i++)//&& (start_x + i < size_x)  !!!
    y[i] = x[start_x + i];

};

void linalg::introduce_segment(const std::vector< std::complex<double> >&x, std::vector< std::complex<double> >& y,  size_t start_y ){//size_y >> size_y


#pragma omp parallel for 
  for(size_t i = 0; i< x.size(); i++)
    y[start_y + i] = x[i];

};  


void linalg::orthogonalize(SparseMatrixType& S, const std::vector< std::complex<double> >& original, std::vector< std::complex<double> >& orthogonalized){


  Eigen::Map<const Eigen::Vector<std::complex<double>, -1>>
    eig_original(original.data(), original.size());
  Eigen::Map<Eigen::Vector<std::complex<double>, -1>>
    eig_orthogonalized(orthogonalized.data(), original.size());




    
  Eigen::ConjugateGradient<  Eigen::SparseMatrix<std::complex<double>,  Eigen::RowMajor, indexType>, Eigen::Lower,  Eigen::DiagonalPreconditioner< std::complex<double>> > solver;
  solver.setTolerance(0.00001); 
  solver.setMaxIterations(20000); 
  solver.compute(S.eigen_matrix());
  eig_orthogonalized = solver.solve(eig_original);
  
  
  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "  max#iterations:" << solver.maxIterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;
  std::cout << "  tolerance :    " << solver.tolerance()      << std::endl;
  std::cout<<  "Vector norm :    " <<eig_original.norm()<<std::endl<<std::endl;    
  
};
