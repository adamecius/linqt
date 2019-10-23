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
#include "mkl.h"

//#include "cusparse_matrices.h"
//#include "cublas_algebra.h"
//#include "cusparse_multiplication.h"


#endif

namespace blas{
    namespace VectorTransform{	
		void axpy(const int Dim,const float *alpha,const float *x, float *y){
			saxpy(Dim,alpha,x,1,y,1);	
			}
		void axpy(const int Dim,const double *alpha,const double *x, double *y){
			daxpy(Dim,alpha,x,1,y,1);	
			}
		void axpy(const int Dim,const std::complex<float> *alpha,const std::complex<float> *x,std::complex<float> *y){
			caxpy(Dim,alpha,x,1,y,1);	
			}
		void axpy(const int Dim,const std::complex<double> *alpha,const std::complex<double> *x, std::complex<double> *y){
			zaxpy(Dim,alpha,x,1,y,1);	
			}
		void dot(const int Dim,const float *x,const float *y,float *result ){
			sdot(Dim,x,1,y,1,result);	
			}
		void dot(const int Dim,const double *x,const double *y,double *result ){
			ddot(Dim,x,1,y,1,result);	
			}
		void dot(const int Dim,const std::complex<float> *x,const std::complex<float> *y,std::complex<float> *result){
			cdotc(Dim,x,1,y,1,result);	
			}
		void dot(const std::complex<double> *x,const std::complex<double> *y, std::complex<double> *result){
			zdotc(Dim,x,1,y,1,result);	
			}
		void swap(const int Dim, float *x, float *y){
			sswap(Dim,x,1,y,1);	
			}
		void swap(const int Dim, double *x, double *y){
			dswap(Dim,x,1,y,1);	
			}
		void swap(const int Dim, std::complex<float> *x, cusp::complex<float> *y){
			cswap(Dim,x,1,y,1);	
			}
		void swap(const int Dim, std::complex<double> *x, cusp::complex<double> *y){
			zswap(Dim,x,1,y,1);	
			}
		void scale( int Dim,const float *alpha,float *x){
			sscal(Dim,alpha,x,1);
			}
		void scale( int Dim,const double *alpha,double *x){
			dscal(Dim,alpha,x,1);
			}
		void scale( int Dim,const float *alpha, cusp::complex<float> *x){
			csscal(Dim,alpha,x,1);
			}
		void scale( int Dim,const std::complex<float> *alpha,cusp::complex<float> *x){
			cscal(Dim,alpha,x,1);
			}
		void scale( int Dim,const std::complex<double>alpha,cusp::complex<double> *x){
			zscal(Dim,alpha,x,1);
			}
		void scale( int Dim,const double *alpha,cusp::complex<double> *x){
			zdscal(Dim,alpha,x,1);
			}

		void nrm2(const int Dim,const float *x,float *result){
			snrm2(Dim,x,1,result);	
			}
		void nrm2(const int Dim,const double *x,double *result){
			dnrm2(Dim,x,1,result);	
			}
		void nrm2(const int Dim,const std::complex<float>  *x,float *result){
			scnrm2(Dim,x,1,result);	
			}
		void nrm2(const int Dim,const std::complex<double> *x,double *result){
			dznrm2(Dim,x,1,result);	
			}
			/***************FUNCION <x|x>****************/
		template <typename Float,typename Scalar>    void normalize( int Dim, Scalar* x){
			Float norm;
			nrm2(handle,Dim,x,&norm);
			norm=1/norm;
			scale(handle,Dim,&norm,x);
			}
		}
	}

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




