#ifndef KPM_HPP
#define KPM_HPP

#include "types_definitions.hpp"
#include "tb_operator.hpp"


void cheb_evolve(	const int n,
					const TBOperator Ham,
					std::vector<kpm::complex> *jn0,
					std::vector<kpm::complex> *jn1
				) 
{
	switch (n)
	{
		case 0:
		{
			Ham.Multiply( 1. , *jn0 , 0 , *jn1 ); 
		}
		break;
		default:
		{
			Ham.Multiply( 2. , *jn0 ,-1 , *jn1 );			

			std::vector<kpm::complex> *jnt = jn1 ;
			jn1 = jn0;
			jn0 = jnt; 
		}
	}  
}


#endif
