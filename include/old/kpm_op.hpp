#ifndef KPM_HPP
#define KPM_HPP

#include "types_definitions.hpp"
#include "kpm_linalg.hpp"


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
