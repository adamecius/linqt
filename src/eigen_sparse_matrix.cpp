#include "sparse_matrix.hpp"


const int MKL_SPMVMUL = 1000;

void SparseMatrixType::ConvertFromCOO(vector<indexType> &rows, vector<indexType> &cols, vector<complex<double> > &vals)
{
	rows_ = vector<indexType>(rows);
	cols_ = vector<indexType>(cols);
	vals_ = vector<complex<double> >(vals);

	std::size_t NNZ = vals_.size();
	

        std::vector<Eigen::Triplet<std::complex<double>,indexType>> triplets;
        triplets.reserve(NNZ);


        for (std::size_t i = 0; i < NNZ; ++i) {
	  triplets.push_back(Eigen::Triplet<std::complex<double>, indexType>(rows_[i], cols_[i], vals_[i]));
        }


        int num_rows = *std::max_element(rows_.begin(), rows_.end()) + 1;
        int num_cols = *std::max_element(cols_.begin(), cols_.end()) + 1;



        Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor, indexType> matrix(num_rows, num_cols);
        matrix.setFromTriplets(triplets.begin(), triplets.end());

	
	setDimensions(num_rows, num_cols);
	std::cout<<"ERROR: Incomplete function for constructing Hamiltonian from COO input matrix. I dont want to be generating extra copies of a huge input."<<std::endl;
        std::terminate();
	//new (&matrix_) Eigen::Map<Eigen::SparseMatrix<complex<double>, Eigen::RowMajor> >(rows_.size(), cols_.size(), NNZ, rows_.data(), cols_.data(), vals_.data());
};

void SparseMatrixType::ConvertFromCSR(vector<indexType> &rowIndex, vector<indexType> &cols, vector<complex<double> > &vals)
{
	rows_ = vector<indexType>(rowIndex);
	cols_ = vector<indexType>(cols);
	vals_ = vector<complex<double> >(vals);

        indexType NNZ = vals_.size();

	
	matrix_ = Eigen::Map<Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType> >(rows_.size()-1, rows_.size()-1, NNZ, rows_.data(), cols_.data(), vals_.data());
	
	//Eigen::SelfAdjointView<Eigen::SparseMatrix<double>, Eigen::Upper> symmetric_view(sparse_matrix); ->Maybe saves half the memory? Couldnt see the documentation

	Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType> transp = Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType>( matrix_.adjoint() ); 
	if( ( transp - matrix_ ).norm() != 0 ){
	  std::cout<<"Hamiltonian has to be symmetrized."<<std::endl;
	  //matrix_ =  transp + matrix_; 
       	}

	/*
	std::cout<<Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>(matrix_).block(0,0,10,10)<<std::endl;

	std::cout<<"HEERE 3: "<<matrix_.innerIndexPtr()[0]<<"  "<<matrix_.innerIndexPtr()[1]<<"  "<<matrix_.innerIndexPtr()[2]<<std::endl;
	std::cout<<"HEERE 3: "<<matrix_.outerIndexPtr()[0]<<"  "<<matrix_.outerIndexPtr()[1]<<"  "<<matrix_.outerIndexPtr()[2]<<std::endl;
	std::cout<<"HEERE 3: "<<matrix_.valuePtr()[0]<<"  "<<matrix_.valuePtr()[1]<<"  "<<matrix_.valuePtr()[2]<<std::endl;

	std::complex<double> errVal=0.0;
	double errCol=0.0, errRow=0.0;
	
	for(indexType i=0;i<(long int)NNZ;i++){
	  errVal += matrix_.valuePtr()[i] - vals_[i];
	  errRow += cols_[i] - matrix_.innerIndexPtr()[i];  
	}

	
	for(indexType i = 0; i < rows_.size(); i++ ){
	  errCol += rows_[i] - matrix_.outerIndexPtr()[i];  
	}

	
	std::cout<<"Errors: val-"<<errVal<<" col-"<<errCol<<" row-"<< errRow<<std::endl;
	*/

	
	setDimensions(rows_.size()-1, rows_.size()-1);

}


void SparseMatrixType::Multiply(const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y)
{
        Eigen::Map<const Eigen::Vector<std::complex<double>, -1>>
	  eig_x(x, numRows());

	Eigen::Map<Eigen::Vector<std::complex<double>, -1>>
	  eig_y(y, numRows());


	eig_y = a * matrix_ * eig_x + b * eig_y; 

	return ;
};

void SparseMatrixType::Multiply(const complex<double> a, const vector< complex<double> >& x, const complex<double> b, vector< complex<double> >& y)
{
        Multiply(a, x.data(), b, y.data() );
        return ;
};



void SparseMatrixType::Rescale(const complex<double> a, const complex<double> b)
{
	//Create Identity Matrix
        Eigen::SparseMatrix<complex<double>, Eigen::RowMajor, indexType> bID(numRows(), numCols()); 


	bID.setIdentity();
        bID *= b;
	
	matrix_ = a * matrix_ + bID;


	return ;
}




void SparseMatrixType::BatchMultiply(const int batchSize, const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y)
{
  	std::cout<<"ERROR: SparseMatrixType:BatchMultiply function not implemented in eigen yet."<<std::endl;
	std::terminate();
	//auto status = mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE,a,Matrix, descr,SPARSE_LAYOUT_COLUMN_MAJOR,x, batchSize, numCols(), b, y, numCols() );
  //assert( status == SPARSE_STATUS_SUCCESS);
}

