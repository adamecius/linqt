

#include "kpm_linalg.hpp"
#include <iostream>

#include "omp.h"
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "mkl_types.h"
#include "mkl.h"

void
KPMLinalg::swap( my::scalar* x , my::scalar* y)
{
	my::scalar* tmp=x;
	x= y;
	y= tmp;
}

void
KPMLinalg::copy( my::scalar* x , my::scalar* y)
{
	#pragma omp parallel for
	for(my::integer i=0;i<eff_dim;i+=shift)
	{
		y[i]=x[i];
	}
}

void
KPMLinalg::scale( 	const my::real alpha,
					my::scalar* x )
{
	#pragma omp parallel  for
	for(my::integer i=0;i<eff_dim;i+=shift)
	{
		x[i]=alpha* x[i] ;
	}

};


void
KPMLinalg::axpy(
				const my::real alpha,
				const my::scalar* x,
				my::scalar* y )
{

	#pragma omp parallel for
	for(my::integer i=0;i<eff_dim;i+=shift)
	{
		y[i]=y[i] + alpha* x[i] ;
	}

};

my::scalar
KPMLinalg::dot(
				const my::scalar* x,
				const my::scalar* y	)
{
	my::scalar sum=0.;

	#pragma omp parallel shared(sum)
	{
		my::scalar  priv_sum = 0.;
		#pragma omp for
		for (my::integer i=0; i<eff_dim;i+=shift )
		{
			priv_sum += std::conj(x[i])*y[i];
		}

		#pragma omp critical
		{
			sum += priv_sum;
		}
	}

	return sum;
};



my::real
KPMLinalg::nrm2( const my::scalar* x)
{
	my::real norm=0.;

	#pragma omp parallel  shared(norm)
	{
		my::real  priv_norm= 0.;
		#pragma omp for
		for (my::integer i=0; i<eff_dim;i+=shift )
		{
			priv_norm += std::norm(x[i]); //returns Re^2+Im^2
		}

		#pragma omp critical
		{
			norm += priv_norm ;
		}
	}
	return sqrt(norm);
};

void
KPMLinalg::normalize( my::scalar* x )
{
	const my::real norm=1./nrm2(x);
	scale(norm,x);
};



