#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX

#include <assert.h>     /* assert */
#include <iostream>
#include <vector>
#include <complex>
using namespace std;
 
class SparseMatrixType
{
  public:
    virtual string matrixType()= 0;
    virtual void Multiply(const complex<double>* a,const complex<double>* x,const complex<double>*  b, complex<double>* y) = 0;
    virtual void Optimize() = 0;
    virtual void ConvertFromCOO(vector<int> &rows,vector<int> &cols,vector<complex<double> > &vals) = 0;
    int numRows(){ return numRows_;};
    int numCols(){ return numCols_;};
    int rank() { return ((this->numRows()>this->numCols())? this->numCols(): this->numRows() ) ; };
    void setDimensions(const int numRows,const int numCols){ numRows_=numRows; numCols_=numCols; };

  protected:
    int numRows_, numCols_;   
};

// MKL LIBRARIES
#define MKL_Complex16 complex<double>
#include "mkl.h"
#include "mkl_spblas.h"
class MKL_SparseType : public SparseMatrixType
{
  public:
  	string matrixType();
	void Multiply(const complex<double>* a,const complex<double>* x,const complex<double>*  b, complex<double>* y);
	void Optimize();
	void ConvertFromCOO(vector<int> &rows,vector<int> &cols,vector<complex<double> > &vals);

	private:
		struct matrix_descr descr;
		sparse_matrix_t Matrix;

};

class SparseMatrixBuilder
{
  public:
    void setSparseMatrix(SparseMatrixType *b)
    {
        _matrix_type = b;
    };

  public:
    void ConstructFromCOO( const int numRows, const int numCols,
                            vector<int>& rows, 
                            vector<int>& cols,
                            vector< complex<double> >& vals)
    {
        std::cout<<"Building: "<<_matrix_type->matrixType()<<std::endl;
        _matrix_type->setDimensions(numRows,numCols);
        _matrix_type->ConvertFromCOO(rows,cols,vals);
        _matrix_type->Optimize();
    };

    SparseMatrixType* _matrix_type;
};




#endif


