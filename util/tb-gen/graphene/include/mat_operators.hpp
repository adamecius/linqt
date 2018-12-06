
#ifndef MAT_OPERATOR_HPP
#define MAT_OPERATOR_HPP

#include "lattice_geometry.h"

inline
void
LinearToSquareConversion(const int dim, const Complex* A, std::vector< std::vector<Complex> >& B )
{
	for(int i=0;i<dim; i++)
	for(int j=0;j<dim; j++)
		B[i][j]=A[j*dim+i];
}


inline
void
SecondEqualFirst(const int dim, const Complex* A, Complex* B )
{
	for(int i=0;i<dim*dim; i++)
		B[i]=A[i];
}


inline
void
MatVecMul(const int dim, const Complex* H,const Complex* x, Complex* y )
{
	Complex out[dim];
	for(int i=0;i<dim; i++)
	{
		out[i]=0;
		for(int j=0;j<dim; j++)
			out[i]+= H[i*dim+j]*x[j];
	}		
	for(int i=0;i<dim; i++)
		y[i]=out[i];
}

inline
void
MatMatMul(const int dim,const Complex* A,const Complex* B, Complex* C )
{
	Complex OUT[dim*dim];
	
	for(int i=0;i<dim; i++)
	for(int j=0;j<dim; j++)
	{
		OUT[i*dim + j]=0;
		for(int k=0;k<dim; k++)
			OUT[i*dim + j]+=A[i*dim+k]*B[k*dim+j];
	}
	for(int n=0;n<dim*dim; n++)
		C[n]=OUT[n];
}


inline
Complex
MatElem(const int dim,const Complex* x,const Complex* H,const Complex* y )
{
	Complex z[dim];
	MatVecMul(dim,H,y,z);
	return Dot(dim,x,z);
}


void
MatAntiCom(const int dim,const Complex* A,const Complex* B,Complex* C)
{
	Complex CL[dim*dim], CR[dim*dim];
	MatMatMul(dim,A,B,CL);
	MatMatMul(dim,B,A,CR);
	
	for(int n=0;n<dim*dim;n++)
		C[n]=(CL[n]+CR[n]);
}

void BasisChange(const int dimk, Complex* M, Complex* U)
{
	//The matrix U is assumed to have its unitary vector stored in 
	// column-wise form
	Complex
	Mnew[dimk*dimk];
	for(int eig0=0;eig0<dimk;eig0++)	   
	for(int eig1=0;eig1<dimk;eig1++)
	{
		Mnew[eig0*dimk + eig1]=0.0;
		for(int m=0; m <dimk; m++)	    
		for(int n=0; n <dimk; n++)
//			Mnew[eig0*dimk + eig1]+=  conj(U[ eig0 *dimk  + m  ])*M[ m*dimk + n]* U[ eig1 *dimk  + n  ];
			Mnew[eig0*dimk + eig1]+=  conj(U[ m *dimk  + eig0  ])*M[ m*dimk + n]* U[ n *dimk  + eig1  ];

	}
	for(int n=0;n<dimk*dimk;n++)
		M[n]=Mnew[n];
};



const Complex IdentMat[DIMK*DIMK]={ 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1} ;

#endif
