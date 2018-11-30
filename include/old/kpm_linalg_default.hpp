#ifndef KPM_LINALG_HPP
#define KPM_LINALG_HPP

#include "types_definitions.hpp"
#include <cstring> //memset

#ifdef USEMKL
#include "mkl.h"
#endif

#ifdef USEOPEN
#include <omp.h>
#endif


namespace kpm
{

	void reset(const size_t dim, kpm::complex* x)
	{
		memset(x,0,dim*sizeof(x[0])); //Clear storage vector
	}

	void dot(const size_t dim, const kpm::complex* x,const kpm::complex* y, kpm::complex&result) 
	{
//		result=0;	
//		for (size_t i=0; i< dim;i++ )
//			result+= std::conj(x[i])*y[i];

		const size_t MAX_THREADS=omp_get_max_threads() ;
		kpm::complex sum_arr[MAX_THREADS];
		#pragma omp parallel shared(sum_arr)
		{
			kpm::complex priv_sum = 0.;
			#pragma omp for
			for (size_t i=0; i<dim;i++ )
				priv_sum += std::conj(x[i])*y[i];
				
				const size_t tidx = omp_get_thread_num();
				sum_arr[tidx]=priv_sum;
		}
		result=0;
		for( size_t n=0;n< MAX_THREADS ; n++)
			result+=sum_arr[n];
	};


	void copy(const size_t dim, kpm::complex* dest, kpm::complex* source ) 
	{
		#pragma omp parallel for
			for (size_t i=0; i<dim;i++ )
				dest[i]=source[i] ;
	};

	void axpy(	const size_t dim			, kpm::complex _a,
				const kpm::complex* source	, kpm::complex* dest) 
	{
		const kpm::complex a=_a;
		#pragma omp parallel for firstprivate(a)
		for (size_t i=0; i< dim;i++ )
			dest[i] = dest[i] + a*source[i] ;
	};

	void axpy(	const size_t dim			, kpm::real _a,
				const kpm::complex* source	, kpm::complex* dest) 
	{
		const kpm::real a=_a;
		#pragma omp parallel for firstprivate(a)
		for(size_t i=0; i< dim; i++ )
			dest[i] = dest[i] + a*source[i] ;

	};


	void reset(const size_t dim, kpm::complex* x)
	{
		memset(x,0,dim*sizeof(x[0])); //Clear storage vector
	}

	void dot(const size_t dim, const kpm::complex* x,const kpm::complex* y, kpm::complex&result) 
	{
//		result=0;	
//		for (size_t i=0; i< dim;i++ )
//			result+= std::conj(x[i])*y[i];

		const size_t MAX_THREADS=omp_get_max_threads() ;
		kpm::complex sum_arr[MAX_THREADS];
		#pragma omp parallel shared(sum_arr)
		{
			kpm::complex priv_sum = 0.;
			#pragma omp for
			for (size_t i=0; i<dim;i++ )
				priv_sum += std::conj(x[i])*y[i];
				
				const size_t tidx = omp_get_thread_num();
				sum_arr[tidx]=priv_sum;
		}
		result=0;
		for( size_t n=0;n< MAX_THREADS ; n++)
			result+=sum_arr[n];
	};


	void copy(const size_t dim, kpm::complex* dest, kpm::complex* source ) 
	{
		#pragma omp parallel for
			for (size_t i=0; i<dim;i++ )
				dest[i]=source[i] ;
	};

	void axpy(	const size_t dim			, kpm::complex _a,
				const kpm::complex* source	, kpm::complex* dest) 
	{
		const kpm::complex a=_a;
		#pragma omp parallel for firstprivate(a)
		for (size_t i=0; i< dim;i++ )
			dest[i] = dest[i] + a*source[i] ;
	};

	void axpy(	const size_t dim			, kpm::real _a,
				const kpm::complex* source	, kpm::complex* dest) 
	{
		const kpm::real a=_a;
		#pragma omp parallel for firstprivate(a)
		for(size_t i=0; i< dim; i++ )
			dest[i] = dest[i] + a*source[i] ;

	};

}

#endif
