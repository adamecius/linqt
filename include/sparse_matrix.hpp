#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX

#include <assert.h> /* assert */
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
using namespace std;

namespace Sparse
{
bool OPERATOR_FromCSRFile(const std::string input, int &dim, vector<int> &columns, vector<int> &rowIndex, vector<complex<double> > &values);
};

class SparseMatrixType
{
public:
  virtual string matrixType() const = 0;
  virtual void Multiply(const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y) = 0;
  virtual void Optimize() = 0;
  virtual void ConvertFromCOO(vector<int> &rows, vector<int> &cols, vector<complex<double> > &vals) = 0;
  virtual void ConvertFromCSR(vector<int> &rowIndex, vector<int> &cols, vector<complex<double> > &vals) = 0;
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

private:
  int numRows_, numCols_;
  string id_;
};

// MKL LIBRARIES
#define MKL_Complex16 complex<double>
#include "mkl.h"
#include "mkl_spblas.h"
class MKL_SparseType : public SparseMatrixType
{
public:
  MKL_SparseType()
  {
    descr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
    descr.mode = SPARSE_FILL_MODE_UPPER;
    descr.diag = SPARSE_DIAG_NON_UNIT;
  }

  string matrixType() const { return "CSR Matrix from MKL Library."; };
  void Multiply(const complex<double> a, const complex<double> *x, const complex<double> b, complex<double> *y);
  void Optimize();
  void ConvertFromCOO(vector<int> &rows, vector<int> &cols, vector<complex<double> > &vals);
  void ConvertFromCSR(vector<int> &rowIndex, vector<int> &cols, vector<complex<double> > &vals);

private:
  struct matrix_descr descr;
  sparse_matrix_t Matrix;
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
    std::cout << "OPERATOR SUCCESSFULLY BUILD" << std::endl;
  }
  SparseMatrixType *_matrix_type;
};

#endif
