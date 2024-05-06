#include "sparse_matrix.hpp"

const int MKL_SPMVMUL = 1000;

void SparseMatrixType::ConvertFromCOO(vector<int> &rows, vector<int> &cols, vector<complex<double> > &vals)
{
	rows_ = vector<int>(rows);
	cols_ = vector<int>(cols);
	vals_ = vector<complex<double> >(vals);

	sparse_matrix_t newMatrix;

	auto status = mkl_sparse_z_create_coo(&newMatrix, SPARSE_INDEX_BASE_ZERO, numRows(), numCols(), rows_.size(), &rows_[0], &cols_[0], &vals_[0]); 
	assert( status == SPARSE_STATUS_SUCCESS);
	
	status = mkl_sparse_convert_csr(newMatrix, SPARSE_OPERATION_NON_TRANSPOSE, &Matrix);
	assert( status == SPARSE_STATUS_SUCCESS);

	status = mkl_sparse_destroy(newMatrix);
	assert( status == SPARSE_STATUS_SUCCESS);

	status = mkl_sparse_z_create_coo(&newMatrix, SPARSE_INDEX_BASE_ZERO, numRows(), numCols(), rows_.size(), &rows_[0], &cols_[0], &vals_[0]);
	assert( status == SPARSE_STATUS_SUCCESS);

	status = mkl_sparse_convert_csr(newMatrix, SPARSE_OPERATION_NON_TRANSPOSE, &Matrix);
	assert( status == SPARSE_STATUS_SUCCESS);

	status = mkl_sparse_destroy(newMatrix);
	assert( status == SPARSE_STATUS_SUCCESS);

	status= mkl_sparse_set_mv_hint (Matrix,SPARSE_OPERATION_NON_TRANSPOSE, descr, MKL_SPMVMUL);
	assert(status == SPARSE_STATUS_SUCCESS);

	status= mkl_sparse_optimize (Matrix);
	assert(status == SPARSE_STATUS_SUCCESS);

};

void SparseMatrixType::ConvertFromCSR(vector<int> &rowIndex, vector<int> &cols, vector<complex<double> > &vals)
{
	rows_ = vector<int>(rowIndex);
	cols_ = vector<int>(cols);
	vals_ = vector<complex<double> >(vals);

	auto status = mkl_sparse_z_create_csr(&Matrix, SPARSE_INDEX_BASE_ZERO, numRows(), numCols(), &rows_[0], &rows_[1], &cols_[0], &vals_[0]) ;
	assert(status == SPARSE_STATUS_SUCCESS);
	
	status= mkl_sparse_set_mv_hint (Matrix,SPARSE_OPERATION_NON_TRANSPOSE, descr, MKL_SPMVMUL);
	assert(status == SPARSE_STATUS_SUCCESS);
 	
	status	= mkl_sparse_optimize (Matrix);
	assert(status == SPARSE_STATUS_SUCCESS);

}

void SparseMatrixType::Multiply(const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y)
{
	auto status = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, a, Matrix, descr, x, b, y);
	assert(status == SPARSE_STATUS_SUCCESS);
	return ;
};

void SparseMatrixType::Multiply(const complex<double> a, const vector< complex<double> >& x, const complex<double> b, vector< complex<double> >& y)
{
	assert(x.size() == y.size());
	auto status = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, a, Matrix, descr, &x[0], b, &y[0]);
	assert(status == SPARSE_STATUS_SUCCESS);
	return ;
};



void SparseMatrixType::Rescale(const complex<double> a,const complex<double> b)
{
	//Create Identity Matrix
	sparse_matrix_t bID; 
	std::vector<int> row(numRows()+1,0); for( int i=0; i < numRows()+1; i++) row[i]=i;
	std::vector< complex<double> >  val(numRows(),b);
	//B := b

	auto status = mkl_sparse_z_create_csr(&bID, SPARSE_INDEX_BASE_ZERO, numRows(), numCols(), &row[0], &row[1], &row[0], &val[0]) ;
	assert(status == SPARSE_STATUS_SUCCESS);

	//COPY the new matrix
	sparse_matrix_t A; 
	status = mkl_sparse_copy (Matrix,descr,&A);
	assert(status == SPARSE_STATUS_SUCCESS);

	//C := alpha*op(A) + B
	//C := a*op(A) + B

	status = mkl_sparse_z_add (SPARSE_OPERATION_NON_TRANSPOSE,A,a,bID, &Matrix);
	assert(status == SPARSE_STATUS_SUCCESS);


	return ;
}




void SparseMatrixType::BatchMultiply(const int batchSize, const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y)
{
	auto status = mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE,a,Matrix, descr,SPARSE_LAYOUT_COLUMN_MAJOR,x, batchSize, numCols(), b, y, numCols() );
	assert( status == SPARSE_STATUS_SUCCESS);
}

