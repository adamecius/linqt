#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX

#include <assert.h> /* assert */
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>

#include<eigen-3.4.0/Eigen/Core>
#include<eigen-3.4.0/Eigen/Sparse>

using std::complex;
using std::vector;
using std::string;

namespace Sparse
{
bool OPERATOR_FromCSRFile(const std::string input, int &dim, vector<int> &columns, vector<int> &rowIndex, vector<complex<double> > &values);
};

class SparseMatrixType_BASE
{
	
	
public:
	typedef complex<double> value_t;
	typedef vector< value_t > vector_t;
	
	  int numRows() { return numRows_; };
	  int numCols() { return numCols_; };
	  int rank() { return ((this->numRows() > this->numCols()) ? this->numCols() : this->numRows()); };
	  void setDimensions(const int numRows, const int numCols)
	  {
		numRows_ = numRows;
		numCols_ = numCols;
	  };
	  void SetID(string id) { id_ = id; }
	  string ID() const { return id_; }
		
	  bool isIdentity(){ return (bool)( ID()=="1"); };

private:
  int numRows_, numCols_;
  string id_;
};


// MKL LIBRARIES
class SparseMatrixType  : public SparseMatrixType_BASE
{
public:
  SparseMatrixType()
  {}
  ~SparseMatrixType() { delete matrix_; }
    
  
  string matrixType() const { return "CSR Matrix from Eigen."; };
  void Multiply(const value_t a, const value_t *x, const value_t b, value_t *y);
  void Multiply(const value_t a, const vector_t& x, const value_t b, vector_t& y);
  void Rescale(const value_t a,const value_t b);
  inline 
  void Multiply(const value_t *x, value_t *y){ Multiply(value_t(1,0),x,value_t(0,0),y);};
  inline 
  void Multiply(const vector_t& x, vector_t& y){ Multiply(value_t(1,0),x,value_t(0,0),y);};


  void BatchMultiply(const int batchSize, const value_t a, const value_t *x, const value_t b, value_t *y);

  void ConvertFromCOO(vector<int> &rows, vector<int> &cols, vector<complex<double> > &vals);
  void ConvertFromCSR(vector<int> &rowIndex, vector<int> &cols, vector<complex<double> > &vals);

  vector<int>* rows() {return &rows_;};
  vector<int>* cols() {return &cols_;};
  vector<complex<double> >* vals(){return &vals_;};

private:
  Eigen::Map<Eigen::SparseMatrix<complex<double>, Eigen::RowMajor> >* matrix_;
  vector<int> rows_;
  vector<int> cols_;
  vector<complex<double> > vals_;
};



class SparseMatrixBuilder
{
public:
  void setSparseMatrix(SparseMatrixType *b)
  {
    _matrix_type = b;
  };

public:
  void BuildOPFromCSRFile(const std::string input)
  {
    vector<int> columns, rowIndex;
    vector<complex<double> > values;
    int dim;
    Sparse::OPERATOR_FromCSRFile(input, dim, columns, rowIndex, values);
    _matrix_type->setDimensions(dim, dim);
    _matrix_type->ConvertFromCSR(rowIndex, columns, values);
    std::cout << "OPERATOR SUCCESSFULLY BUILT" << std::endl;
  }
  SparseMatrixType *_matrix_type;
};

#endif
