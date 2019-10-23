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
#ifndef UTILIDADES_H
#define UTILIDADES_H
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
//#include "cusparse_matrices.h"
//#include "cublas_algebra.h"
//#include "cusparse_multiplication.h"


#endif



namespace cusp {
    namespace VectorTransform{
        
        template <typename T> struct saxmy_functor{
            T a;
            saxmy_functor(T a) : a(a) {}
            __host__ __device__
            T operator()(T x, T y) const{return a* x-y;}
        };
        template <typename Vector>  void saxmy(const typename Vector::value_type a, Vector& x, const Vector& y){
            typedef typename Vector::value_type T;
            thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), saxmy_functor<T>(a));
        }
        
        /***************FUNCION Y <- A * Y +X****************/
        template <typename T> struct saxpy_functor{
            T a;
            saxpy_functor(T a) : a(a) {}
            __host__ __device__
            T operator()(T x, T y) const{return a* x+y;}
        };
        template <typename Vector>  void saxpy(const typename Vector::value_type a,  Vector& x,const Vector& y){
            typedef typename Vector::value_type T;
            thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), saxpy_functor<T>(a));
        }

        /***************FUNCION Y <-  Y +A*X****************/
        template <typename T> struct saypx_functor{
            T a;
            saypx_functor(T a) : a(a) {}
            __host__ __device__
            T operator()(T x, T y) const{return  a*x+y;}
        };
        template <typename Vector>  void saypx(const typename Vector::value_type a,const Vector& x, Vector& y){
            typedef typename Vector::value_type T;
            thrust::transform(x.begin(), x.end(), y.begin(), y.begin(), saypx_functor<T>(a));
        }
        
        /***************FUNCION <X,Y>****************/
        template <typename T> struct dot_functor{
            __host__ __device__
            T operator()(T x, T y) const{return conj(x)*y;}
        };
        template <typename Vector>  const typename Vector::value_type dot(const Vector& x,const  Vector& y){
            Vector z;
            z=x;
            typedef typename Vector::value_type T;
            thrust::transform(x.begin(), x.end(),y.begin(), z.begin(), dot_functor<T>());
            return thrust::reduce(z.begin(), z.end());
        }
        
        /***************FUNCION a|x>****************/
        template <typename T> struct scale_functor{
            T a;
            scale_functor(T a) : a(a) {}
            __host__ __device__
            T operator()(T x) const{return a*x;}
        };
        template <typename Matriz> void matrixscale(const typename Matriz::value_type a, Matriz& M){
            // M <- a * M
            typedef typename Matriz::value_type T;
            thrust::transform(M.values.begin(), M.values.end(), M.values.begin(), scale_functor<T>(a));
        }

        template <typename Vector> void scale(const typename Vector::value_type a, Vector& V){
            // M <- a * M
            typedef typename Vector::value_type T;
            thrust::transform(V.begin(), V.end(), V.begin(), scale_functor<T>(a));
        }
        
        /***************FUNCION <x|x>****************/
        template <typename T> struct norm2_functor{
            __host__ __device__
            T operator()(T x) const{return conj(x)*x;}
        };
        template <typename Vector>  const typename Vector::value_type norm2(Vector& x){
            typedef typename Vector::value_type T;
            T norm;
            Vector z;
            z=x;
            thrust::transform(x.begin(), x.end(), z.begin(), norm2_functor<T>());
            return cusp::sqrt(thrust::reduce(z.begin(), z.end()));
        }
        /***************FUNCION <x|x>****************/
        template <typename Vector>    void normalize(Vector& x){
            Vector z;
            typedef typename Vector::value_type T;
            T norm;
            T ONE=1;
            z=x;
            thrust::transform(x.begin(), x.end(), z.begin(), norm2_functor<T>());
            norm=thrust::reduce(z.begin(), z.end());
            norm=ONE/cusp::sqrt(norm);
            thrust::transform(x.begin(), x.end(), x.begin(), scale_functor<T>(norm));
        }
        
        
        
    }
    
}

//This function will show how much left to finish one or many loops
static inline void	cycletime(int ncycle){
    

	static	time_t start,end;											//Static variables that will save the time
	static	double total=0;												//Total time accumulated	
	static	int   tint  =1;												//time in second to show the results					
	static	int ncalls	  ;													//number of cycles
	static	int i=0;														//
    
    if(ncycle==-1) {i=0; ncalls=0; total=0; tint  =1; }else{							//This will allow to reset the variables if necessary
	if(i==0){start=clock(); i=1; }										//This line will catch the time at the begining of the cycle (i=0)
	else {																       
        end=clock();   i=0;												//this will catch the time at the end of cycle
        total =total+(double)((double)(end-start)/(double)CLOCKS_PER_SEC);	//this will calculate the diference  in seconds
        ncalls=ncalls+1;												//this will store the  number of cycles
		//std::cout<<"time enlapsep="<<total<<" ncalls="<<ncalls<<" "<<ncycle<<std::endl;
        if(total>=tint){
			tint=tint+1;
			double tcycle=total/(double)ncalls; 						//time per cycle	
			double rest  = tcycle*((double)(ncycle-ncalls));			//time until end
			if(rest <60 ){												//seconds
				printf ("The cycle mean time is  %.2f ms. The whole set wil end in %.2f sec, %d%% completed", 1000*tcycle,rest,(int)(ncalls*100.0f/ncycle));
				printf("\n\033[F\033[J");
				}
			if(rest >=60&&rest <3600){									//minutes
				printf ("The cycle mean time is  %.2f s. The whole set wil end in %.2f min, %d%% completed", tcycle,rest/60.0f,(int)(ncalls*100.0f/ncycle));
				printf("\n\033[F\033[J");
				}
			if(rest >=3600&&rest <86400){								//hours
				printf ("The cycle mean time is  %.2f s. The whole set wil end in %.2f hours, %d%% completed", tcycle,rest/3600.0f,(int)(ncalls*100.0f/ncycle));
				printf("\n\033[F\033[J");
				}
			if(rest >=86400				){								//days
				printf ("The cycle mean time is  %.2f min. The whole set wil end in %.2f days, %d%% completed", tcycle/60.0f,rest/86400.0f,(int)(ncalls*100.0f/ncycle));
				printf("\n\033[F\033[J");
				}

            }
	}
							}
}

