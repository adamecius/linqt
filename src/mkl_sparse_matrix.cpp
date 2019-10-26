#include "sparse_matrix.hpp"

void MKL_SparseType::ConvertFromCOO(vector<int> &rows, vector<int> &cols, vector<complex<double>> &vals)
{
	rows_ = vector<int>(rows);
	cols_ = vector<int>(cols);
	vals_ = vector<complex<double>>(vals);

	sparse_matrix_t newMatrix;
	assert(mkl_sparse_z_create_coo(&newMatrix, SPARSE_INDEX_BASE_ZERO, numRows(), numCols(), rows_.size(), &rows_[0], &cols_[0], &vals_[0]) == SPARSE_STATUS_SUCCESS);
	assert(mkl_sparse_convert_csr(newMatrix, SPARSE_OPERATION_NON_TRANSPOSE, &Matrix) == SPARSE_STATUS_SUCCESS);
	assert(mkl_sparse_destroy(newMatrix) == SPARSE_STATUS_SUCCESS);
};

void MKL_SparseType::ConvertFromCSR(vector<int> &rowIndex, vector<int> &cols, vector<complex<double>> &vals)
{
	rows_ = vector<int>(rowIndex);
	cols_ = vector<int>(cols);
	vals_ = vector<complex<double>>(vals);

	assert(mkl_sparse_z_create_csr(&Matrix, SPARSE_INDEX_BASE_ZERO, numRows(), numCols(), &rows_[0], &rows_[1], &cols_[0], &vals_[0]) == SPARSE_STATUS_SUCCESS);
}

void MKL_SparseType::Multiply(const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y)
{
	for (int i = 0; i < 10; i++)
		std::cout << rows_[i] << " " << vals_[i] << std::endl;

	std::cout << "I am multiplying here" << a << " " << b << " " << this->ID() << std::endl;
	assert(mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, a, Matrix, descr, x, b, y) == SPARSE_STATUS_SUCCESS);
	std::cout << "Finished multiplying" << a << " " << b << std::endl;
};

void MKL_SparseType::Optimize(){

	//	std::cout<<*this->rows_[1]<<" "<<*this->rows_[2]<<std::endl;
	/*		mkl_sparse_set_mv_hint(A,SPARSE_OPERATION_NON_TRANSPOSE,descr, (MKL_INT)numMul );
		mkl_sparse_set_dotmv_hint(A,SPARSE_OPERATION_NON_TRANSPOSE,descr, (MKL_INT)numMul );
		mkl_sparse_optimize (A);*/

};
