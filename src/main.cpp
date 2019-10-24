// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include "sparse_matrix.hpp"



int main(int argc, char *argv[])
{
	const int numRows = 10;
	const int numCols = numRows;
	const int nnz = 2*numRows;
	std::vector<int> rows(nnz,0); 
	std::vector<int> cols(nnz,0);
	std::vector<complex<double> > vals(nnz);

	for(int i =0; i < numRows;i++ )
	{rows[2*i]=i; rows[2*i+1]=i; }

	for(int i =0; i < numRows;i++ )
	{cols[2*i]=0 ; cols[2*i+1]=1; }

	for(int i =0; i < nnz;i++ )
		vals[i] = rand()%3;


	MKL_SparseType H;
	SparseMatrixBuilder builder;
	builder.setSparseMatrix(&H);
	builder.ConstructFromCOO(numRows,numCols,rows,cols,vals);

	vector < complex<double> > vecIn(numRows); 	for(int i = 0 ; i < numRows; i++ ) vecIn[i] =10;
	vector < complex<double> > vecOut(numRows);   for(int i = 0 ; i < numRows; i++ ) vecOut[i]=0;
	const complex<double> a=10.0; const complex<double> b=0.0;
	H.Multiply(&a, &vecIn[0], &b, &vecOut[0] );
	for(int i = 0 ; i < numRows; i++ ) std::cout<<vecOut[i]<<std::endl;	
	return 0;
}
