//
//  utilidades.h
//
//
//  Created by Jose Hugo Garcia Aguilar on 14/09/12.
//
//
// Este es un archivo fuente que espero acumule progresivamente funciones
// de uso general para varios programas en C cientificos.

#pragma once
#ifndef LINEAR_CUDA_UTILITIES_H
#define LINEAR_CUDA_UTILITIES_H
//Standar C/C++ headers//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <iostream>
#include <cusp/detail/random.h>
//Cusp headers//
#include <cusp/io/matrix_market.h>
#include <cusp/multiply.h>
#include <cusp/print.h>
#include <cusp/array2d.h>
#include <cusp/coo_matrix.h>
#include <cusp/csr_matrix.h>
#include <cusp/dia_matrix.h>
#include <cusp/ell_matrix.h>
#include <cusp/hyb_matrix.h>
#include <cusp/elementwise.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/device_vector.h>
//Extra Headers//
#include <cusp/complex.h>
#include <cusp/cmath.h>
#include <thrust/reduce.h>
#include <vector>
#include "cublas_v2.h"
#include "cusparse_v2.h"

//#include "cusparse_matrices.h"
//#include "cublas_algebra.h"
//#include "cusparse_multiplication.h"


#endif



namespace cublas {
    namespace VectorTransform{
	void axpy(cublasHandle_t& handle,const int Dim,const float *alpha,const float *x, float *y){
		cublasSaxpy(handle,Dim,alpha,x,1,y,1);	
		}
	void axpy(cublasHandle_t& handle,const int Dim,const double *alpha,const double *x, double *y){
		cublasDaxpy(handle,Dim,alpha,x,1,y,1);	
		}
	void axpy(cublasHandle_t& handle,const int Dim,const cusp::complex<float> *alpha,const cusp::complex<float> *x,cusp::complex<float> *y){
		cublasCaxpy(handle,Dim,alpha,x,1,y,1);	
		}
	void axpy(cublasHandle_t& handle,const int Dim,const cusp::complex<double> *alpha,const cusp::complex<double> *x, cusp::complex<double> *y){
		cublasZaxpy(handle,Dim,alpha,x,1,y,1);	
		}
	void dot(cublasHandle_t& handle,const int Dim,const float *x,const float *y,float *result ){
		cublasSdot(handle,Dim,x,1,y,1,result);	
		}
	void dot(cublasHandle_t& handle,const int Dim,const double *x,const double *y,double *result ){
		cublasDdot(handle,Dim,x,1,y,1,result);	
		}
	void dot(cublasHandle_t& handle,const int Dim,const cusp::complex<float> *x,const cusp::complex<float> *y,cusp::complex<float> *result){
		cublasCdotc(handle,Dim,x,1,y,1,result);	
		}
	void dot(cublasHandle_t& handle,const int Dim,const cusp::complex<double> *x,const cusp::complex<double> *y, cusp::complex<double> *result){
		cublasZdotc(handle,Dim,x,1,y,1,result);	
		}
	void swap(cublasHandle_t& handle,const int Dim, float *x, float *y){
		cublasSswap(handle,Dim,x,1,y,1);	
		}
	void swap(cublasHandle_t& handle,const int Dim, double *x, double *y){
		cublasDswap(handle,Dim,x,1,y,1);	
		}
	void swap(cublasHandle_t& handle,const int Dim, cusp::complex<float> *x, cusp::complex<float> *y){
		cublasCswap(handle,Dim,x,1,y,1);	
		}
	void swap(cublasHandle_t& handle,const int Dim, cusp::complex<double> *x, cusp::complex<double> *y){
		cublasZswap(handle,Dim,x,1,y,1);	
		}
	void scale(cublasHandle_t& handle, int Dim,const float *alpha,float *x){
		cublasSscal(handle,Dim,alpha,x,1);
		}
	void scale(cublasHandle_t& handle, int Dim,const double *alpha,double *x){
		cublasDscal(handle,Dim,alpha,x,1);
		}
	void scale(cublasHandle_t& handle, int Dim,const float *alpha, cusp::complex<float> *x){
		cublasCsscal(handle,Dim,alpha,x,1);
		}
	void scale(cublasHandle_t& handle, int Dim,const cusp::complex<float> *alpha,cusp::complex<float> *x){
		cublasCscal(handle,Dim,alpha,x,1);
		}
	void scale(cublasHandle_t& handle, int Dim,const cusp::complex<double> *alpha,cusp::complex<double> *x){
		cublasZscal(handle,Dim,alpha,x,1);
		}
	void scale(cublasHandle_t& handle, int Dim,const double *alpha,cusp::complex<double> *x){
		cublasZdscal(handle,Dim,alpha,x,1);
		}

	void nrm2(cublasHandle_t& handle,const int Dim,const float *x,float *result){
		cublasSnrm2(handle,Dim,x,1,result);	
		}
	void nrm2(cublasHandle_t& handle,const int Dim,const double *x,double *result){
		cublasDnrm2(handle,Dim,x,1,result);	
		}
	void nrm2(cublasHandle_t& handle,const int Dim,const cusp::complex<float>  *x,float *result){
		cublasScnrm2(handle,Dim,x,1,result);	
		}
	void nrm2(cublasHandle_t& handle,const int Dim,const cusp::complex<double> *x,double *result){
		cublasDznrm2(handle,Dim,x,1,result);	
		}
        /***************FUNCION <x|x>****************/
	template <typename Float,typename Scalar>    void normalize(cublasHandle_t& handle, int Dim, Scalar* x){
		Float norm;
		nrm2(handle,Dim,x,&norm);
		norm=1/norm;
		scale(handle,Dim,&norm,x);
		}
	}
}

namespace cusparse{
		void convertX2Hyb(cusp::coo_matrix<int,float,cusp::host_memory>& H,int DD,cusparseHandle_t& cusparse_handle,cusparseMatDescr_t&	cusparse_descr,cusparseHybMat_t&	MatHyb){
		typedef float	Scalar;
		typedef int		Integer;
		cusp::csr_matrix<Integer,Scalar,cusp::device_memory>		HCSR=H;
		Integer *ptr_csrrow=thrust::raw_pointer_cast(HCSR.row_offsets.data());
		Integer *ptr_csrcol=thrust::raw_pointer_cast(HCSR.column_indices.data());
		Scalar  *ptr_csrval=thrust::raw_pointer_cast(HCSR.values.data());
		cusparseScsr2hyb(cusparse_handle,DD,DD,cusparse_descr,ptr_csrval,ptr_csrrow,ptr_csrcol,MatHyb,0,CUSPARSE_HYB_PARTITION_AUTO);
		}
		void convertX2Hyb(cusp::coo_matrix<int,double,cusp::host_memory>& H,int DD,cusparseHandle_t& cusparse_handle,cusparseMatDescr_t&	cusparse_descr,cusparseHybMat_t&	MatHyb){
		typedef double	Scalar;
		typedef int		Integer;
		cusp::csr_matrix<Integer,Scalar,cusp::device_memory>		HCSR=H;
		Integer *ptr_csrrow=thrust::raw_pointer_cast(HCSR.row_offsets.data());
		Integer *ptr_csrcol=thrust::raw_pointer_cast(HCSR.column_indices.data());
		Scalar  *ptr_csrval=thrust::raw_pointer_cast(HCSR.values.data());
		cusparseDcsr2hyb(cusparse_handle,DD,DD,cusparse_descr,ptr_csrval,ptr_csrrow,ptr_csrcol,MatHyb,0,CUSPARSE_HYB_PARTITION_AUTO);
		}
		void convertX2Hyb(cusp::coo_matrix<int,cusp::complex<float>,cusp::host_memory>& H,int DD,cusparseHandle_t& cusparse_handle,cusparseMatDescr_t&	cusparse_descr,cusparseHybMat_t&	MatHyb){
		typedef cusp::complex<float> Scalar;
		typedef int Integer;
		cusp::csr_matrix<Integer,Scalar,cusp::device_memory>		HCSR=H;
		Integer *ptr_csrrow=thrust::raw_pointer_cast(HCSR.row_offsets.data());
		Integer *ptr_csrcol=thrust::raw_pointer_cast(HCSR.column_indices.data());
		Scalar  *ptr_csrval=thrust::raw_pointer_cast(HCSR.values.data());
		cusparseCcsr2hyb(cusparse_handle,DD,DD,cusparse_descr,ptr_csrval,ptr_csrrow,ptr_csrcol,MatHyb,0,CUSPARSE_HYB_PARTITION_AUTO);
		}
		void convertX2Hyb(cusp::coo_matrix<int,cusp::complex<double>,cusp::host_memory>& H,int DD,cusparseHandle_t& cusparse_handle,cusparseMatDescr_t&	cusparse_descr,cusparseHybMat_t&	MatHyb){
		typedef cusp::complex<double> Scalar;
		typedef int Integer;
		cusp::csr_matrix<Integer,Scalar,cusp::device_memory>		HCSR=H;
		Integer *ptr_csrrow=thrust::raw_pointer_cast(HCSR.row_offsets.data());
		Integer *ptr_csrcol=thrust::raw_pointer_cast(HCSR.column_indices.data());
		Scalar  *ptr_csrval=thrust::raw_pointer_cast(HCSR.values.data());
		cusparseZcsr2hyb(cusparse_handle,DD,DD,cusparse_descr,ptr_csrval,ptr_csrrow,ptr_csrcol,MatHyb,0,CUSPARSE_HYB_PARTITION_AUTO);
		}
		void Multiply(cusparseHandle_t& cusparse_handle,const cusparseMatDescr_t& cusparse_descr,const float *alpha					,const cusparseHybMat_t&	MatHyb,float *x, float* beta, float* y){
			cusparseShybmv(cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,alpha,cusparse_descr,MatHyb,x,beta,y);
               		}
		void Multiply(cusparseHandle_t& cusparse_handle,const cusparseMatDescr_t& cusparse_descr,const double *alpha				,const cusparseHybMat_t&	MatHyb,double *x, double* beta, double* y){
			cusparseDhybmv(cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,alpha,cusparse_descr,MatHyb,x,beta,y);
               		}
		void Multiply(cusparseHandle_t& cusparse_handle,const cusparseMatDescr_t& cusparse_descr,const cusp::complex<float> *alpha	,const cusparseHybMat_t&	MatHyb,cusp::complex<float> *x, cusp::complex<float>* beta, cusp::complex<float>* y){
			cusparseChybmv(cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,alpha,cusparse_descr,MatHyb,x,beta,y);
               		}
		void Multiply(cusparseHandle_t& cusparse_handle,const cusparseMatDescr_t& cusparse_descr,const cusp::complex<double> *alpha	,const cusparseHybMat_t&	MatHyb,cusp::complex<double> *x, cusp::complex<double>* beta, cusp::complex<double>* y){
			cusparseZhybmv(cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,alpha,cusparse_descr,MatHyb,x,beta,y);
               		}

}

namespace linalg{
	
	template <typename Matrix,typename Floating>	void SpectralBounds(Matrix& H, Floating& Emin, Floating& Emax,const Floating Tolerance)
	{
		typedef typename Matrix::value_type Scalar;
		typedef typename Matrix::index_type Integer;
		typedef typename Matrix::memory_space MemorySpace;

		cusparseHandle_t	cusparse_handle;		//Handle of cusparse
		cusparseMatDescr_t	cusparse_descr;			//Cusparse's Matrix Descriptor  
		cusparseHybMat_t	GH;						//Cusparse's Hybrid Matrix
		cublasHandle_t		cublas_handle;			//Handle of cublas
		cusparseCreate			(&cusparse_handle);	//Initialization of the cusparse's handler
		cusparseCreateMatDescr	(&cusparse_descr);	//Initialization of the descriptor			
		cusparseCreateHybMat	(&GH);				//Initialization of the Hybrid Matrix
		cublasCreate			(&cublas_handle);	//Initialization of the cublas's handler  
		Integer DD=H.num_cols;					//Dimensions for the matrix DD*DD
		cusp::coo_matrix<Integer, Scalar, cusp::host_memory> tempH=H;
		cusparse::convertX2Hyb(tempH,DD,cusparse_handle,cusparse_descr,GH);
	
		/*This method is design to find the spectral bounds of any give simmetric matrix
		 * with at least one eigenvalue passing through zero
		 * we will use the power method in order to find the largest eigenvalue in absolute value
		 * and then translate the eigenvalues to find the largest eigenvalue with oposite sign
		 */
		cusp::array1d	<Scalar, cusp::host_memory>   X(DD,0);
		cusp::array1d	<Scalar, cusp::device_memory> GX(DD,0);
		cusp::array1d   <Scalar, cusp::device_memory> GY(DD,0);
		Scalar* pGX=thrust::raw_pointer_cast(GX.data());
		Scalar* pGY=thrust::raw_pointer_cast(GY.data());
		Scalar zero =	 0.0f;
		Scalar one	=	 1.0f;
		Scalar diff	=	-1.0f;
		/***The method parameters***/
		Integer it;				//Iteration Variable
		Integer itmax=100;		//Iteration max
		Integer NRandom=10;		//Number of test vecto we will try
		Floating ConvRate;		//Convergence Rate Variable
		Scalar Etemp=0.0;		//Temporal Variable for the energy
		Scalar E0   =0.0;		//Temporal Variable for the energy
		/***We start the iteration cycle, over the test vectors***/
		for(int r=0;r<NRandom;r++){
			//Select a random vector
			for(int i=0;i<DD;i++) X[i]=(Floating)rand();
			GX =X;
			//Normalize it
			cublas::VectorTransform::normalize<Floating>(cublas_handle,DD,pGX);//Here we normalize j0 and then we pass it to jn
			//Start The convergence cycle//
			ConvRate=1;
			it=0;
			while( (ConvRate>=Tolerance)&& (it<= itmax) ){
				it=it+1;
				//We first compute the vector Y=H X
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pGX,&zero,pGY);
				//Normalize it
				cublas::VectorTransform::normalize<Floating>(cublas_handle,DD,pGY);
				//And then we calculate the quantity X=Y-X
				cublas::VectorTransform::axpy(cublas_handle,DD,&diff,pGY, pGX);
				//The convergence parameter we choose is the norm |Y-X|
				cublas::VectorTransform::nrm2(cublas_handle,DD,pGX,&ConvRate);
				//Notice that because the eigenvalues are defined within a phase
				//it may happen that Y and X be equal up to a phase of -1
				//in this case the best thing to do is just take it out
				ConvRate=abs(ConvRate-2*((int)(ConvRate/(sqrt(2)))));
				GX=GY;
			}
		/**Calculating the First Bound**/
		cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pGX,&zero,pGY);
		cublas::VectorTransform::dot(cublas_handle,DD,pGX,pGY,&Etemp);

		//Here we substitute either by the minimal or maximal bound
		if( Etemp.real()/fabs(Etemp.real()) >0 ) if(Emax < Etemp.real())Emax=Etemp.real();
		if( Etemp.real()/fabs(Etemp.real()) <0 ) if(Emin > Etemp.real())Emin=Etemp.real();
		if(fabs(E0.real()) < fabs(Etemp.real())) E0=Etemp;
		}

		//Now we shift the hamiltonian		
		for(int k=0;k<tempH.num_entries;k++)
			if(tempH.row_indices[k]==tempH.column_indices[k])
				tempH.values[k]=tempH.values[k]-E0;
		cusparse::convertX2Hyb(tempH,DD,cusparse_handle,cusparse_descr,GH);

		//Start the Iterations again
		for(int r=0;r<NRandom;r++){
			//Select a random vector
			for(int i=0;i<DD;i++) X[i]=(Floating)rand();
			GX =X;
			//Normalize it
			cublas::VectorTransform::normalize<Floating>(cublas_handle,DD,pGX);//Here we normalize j0 and then we pass it to jn
			//Start The convergence cycle//
			ConvRate=1;
			it=0;
			while( (ConvRate>=Tolerance)&& (it<= itmax) ){
				it=it+1;
				//We first compute the vector Y=H X
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pGX,&zero,pGY);
				//Normalize it
				cublas::VectorTransform::normalize<Floating>(cublas_handle,DD,pGY);
				//And then we calculate the quantity X=Y-X
				cublas::VectorTransform::axpy(cublas_handle,DD,&diff,pGY, pGX);
				//The convergence parameter we choose is the norm |Y-X|
				cublas::VectorTransform::nrm2(cublas_handle,DD,pGX,&ConvRate);
				//Notice that because the eigenvalues are defined within a phase
				//it may happen that Y and X be equal up to a phase of -1
				//in this case the best thing to do is just take it out
				ConvRate=abs(ConvRate-2*((int)(ConvRate/(sqrt(2)))));
				GX=GY;
			}
		/**Calculating the First Bound**/
		cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pGX,&zero,pGY);
		cublas::VectorTransform::dot(cublas_handle,DD,pGX,pGY,&Etemp);
		Etemp=Etemp;
		//Here we substitute either by the minimal or maximal bound
		if( Etemp.real()/fabs(Etemp.real()) >0 ) if(Emax < Etemp.real())Emax=Etemp.real();
		if( Etemp.real()/fabs(Etemp.real()) <0 ) if(Emin > Etemp.real())Emin=Etemp.real();
		std::cout<<Etemp<<" "<<Emax<<" "<<Emin<<" "<<E0<<std::endl;
		}
		Emin=floorf(Emin)+0.7;
		Emax=ceilf (Emax)-0.7;
		std::cout<<"Emin="<<Emin<<" Emax="<<Emax<<std::endl;
		
		
        
    }
}    
