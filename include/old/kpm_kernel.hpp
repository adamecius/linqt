#ifndef KPM_KERNEL_HPP
#define KPM_KERNEL_HPP

#include "types_definitions.hpp"
#include "structural.hpp"
#include <cassert>


namespace kpm
{
	kpm::real
	JacksonKer( const kpm::real M,const kpm::real m)
	{
		assert( M >0 && m>=0 );
		const kpm::real theta= M_PI/(M+1);
		return ( (M-m+1.)*cos( m*theta)+sin(m*theta)/tan(theta) )/(M+1.);
	}

	kpm::real
	JacksonKerRes( const kpm::real M, const kpm::real x)
	{
		assert( M>=0&&x<1 && x>-1 );
		if( M== 0 )
			return 1;	//When M=0, return the largest possible kernel
		return sqrt( ( M-x*x*(M-1.) )*( 1. - cos(2.*M_PI/(M+1.)) )/(2.*(M+1.)) );
	}

	size_t
	JacksonKerMom4Res( const kpm::real eta, const kpm::real x)
	{
		assert( eta<=1&& eta >0 && x<1 && x>-1 );	

		size_t MaxMom=0; //When eta=1, return the smalles possible Moments 
		if( eta==1 )
			return MaxMom;

		MaxMom=1; 
		kpm::real res=JacksonKerRes( MaxMom,x);
		while ( JacksonKerRes( MaxMom,x) >= eta )
			MaxMom++;
		
		return MaxMom;
	}

	
	kpm::real
	LorentzKer( const kpm::real M,const kpm::real m,const kpm::real lambda=4)
	{
		assert( M >0 && m>=0 && lambda>0);
		const kpm::real theta= 1.0 - (kpm::real)m/(kpm::real)M;
		return sinh( lambda*theta)/sinh(lambda);
	}

	kpm::real
	LorentzKerRes( const kpm::real M, const kpm::real x,const kpm::real lambda=4)
	{
		assert( M >=0 && x<1 && x>-1 && lambda>0);
		if( M== 0 )
			return 1;	//When M=0, return the largest possible kernel

		return  lambda/M;
	}

	size_t
	LorentzKerMom4Res( const kpm::real eta, const kpm::real x,const kpm::real lambda=4)
	{
		assert( eta<=1&& eta >0 && x<1 && x>-1 && lambda>0);
		size_t MaxMom=0; //When eta=1, return the smalles possible Moments 
		if( eta==1 )
			return MaxMom;
		MaxMom = (lambda/eta);
		return MaxMom;
	}




	
	
}
#endif
