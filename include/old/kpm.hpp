#ifndef KPM_HPP
#define KPM_HPP

#include "types_definitions.hpp"
#include "kpm_linalg.hpp"

#ifdef USE_DENSE
#include "dense_tb_operator.hpp"
#else
#include "sparse_tb_operator.hpp"
#endif

void cheb_evolve(	const int n,
					const TBOperator& Ham,
					kpm::complex* &jn0,
					kpm::complex* &jn1
				) 
{
	switch (n)
	{
		case 0:
		{
			Ham.Multiply( 1. , jn0 , 0 , jn1 ); 
			break;
		}

		default:
		{
			kpm::complex *jnt;
			Ham.Multiply( 2. , jn1 ,-1. , jn0 );			
			jnt= jn1 ;
			jn1= jn0;
			jn0= jnt; 
			break;
		}
	}  
}
#endif

/*
#include <gsl/gsl_sf_bessel.h>
void ChebyshevWeightedIteration(	const size_t Dim, 
									const TBOperator& Ham, 
									const size_t Mom,							
									const kpm::complex (*WF)( kpm::real , kpm::real, int),
									const kpm::real param,
									kpm::complex* jm_array,
									kpm::complex* Input)
{
	kpm::complex 
	*jm0=&jm_array[0*Dim],
	*jm1=&jm_array[1*Dim];
	kpm::copy(Dim, jm0 , Input);			//pass init vector to the iteration
	kpm::reset( Dim , Input );		//set the array to zero
	for(int  m=0;m< Mom ; m++ )
	{
		kpm::complex weight = WF( param, m, Mom );
		cheb_evolve( m , Ham, jm0, jm1 );
		kpm::axpy( Dim, weight, jm0, Input );
	}
}
void ChebyshevTimeEvolution(const size_t Dim, 
							const TBOperator& Ham, 
							const size_t Mom,							
							const kpm::real scalFact,  
							kpm::complex* jm_array,
							kpm::complex* Input)
{
	kpm::complex 
	*jm0=&jm_array[0*Dim],
	*jm1=&jm_array[1*Dim];
	kpm::copy(Dim, jm0,Input);			//pass init vector to the iteration
	kpm::reset( Dim , Input);		//set the array to zero
	for(int  m=0;m< Mom ; m++ )
	{
		kpm::complex
		weight = 2.0*gsl_sf_bessel_Jn(m, scalFact) * pow(-I, m );
		if( m ==0 )
			weight*= 0.5;
		cheb_evolve( m , Ham, jm0, jm1 );
		kpm::axpy( Dim, weight, jm0, Input);
	}
}

#endif
*/
