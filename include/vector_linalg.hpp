#ifndef QT_LINALG_HPP
#define QT_LINALG_HPP

#include "types_definitions.hpp"
#include <cstring> //memcpy

namespace qt
{

	inline
	void reset(const qt::dimension dim, qt::complex* x)
	{
		memset(x,0,dim*sizeof(x[0])); //Clear storage vector
	}

	inline
	void copy(const qt::dimension dim, qt::complex* dest, qt::complex* source ) 
	{
		memcpy( dest, source, sizeof( qt::complex)*dim ) ;
	};

	inline
	void dot(const qt::dimension dim, const qt::complex* x,const qt::complex* y, qt::complex& result) 
	{
		result=0.;
		for (qt::index i=0; i< dim;i++ )
		{
			result.real(result.real()+ x[i].real()*y[i].real()+x[i].imag()*y[i].imag() );
			result.imag(result.imag()+ x[i].real()*y[i].imag()-x[i].imag()*y[i].real() );
		}
	};


	inline
	void axpy(	const qt::dimension dim			, qt::complex _a,
				const qt::complex* source	, qt::complex* dest) 
	{
		const qt::complex a=_a;
		for (qt::index i=0; i< dim;i++ )
			dest[i] = dest[i] + a*source[i] ;
	};

	inline
	void axpy(	const qt::dimension dim			, qt::real _a,
				const qt::complex* source	, qt::complex* dest) 
	{
		const qt::real a=_a;
		for(qt::index i=0; i< dim; i++ )
			dest[i] = dest[i] + a*source[i] ;
	};

}

#endif
