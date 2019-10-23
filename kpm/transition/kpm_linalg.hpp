/*
 * kpm_linalg.hpp
 *
 *  Created on: 24/08/2016
 *      Author: jgarcia
 */

#ifndef KPM_LINALG_HPP_
#define KPM_LINALG_HPP_


#include "types_definitions.hpp"
#include "kpm_parallel.hpp"


class KPMLinalg{

	///The Main constructor
public:
	KPMLinalg(const my::integer _dim,const my::integer _shift):
		dim(_dim), shift(_shift), eff_dim(_dim*_shift)
	{}

	my::scalar
	dot( const my::scalar* x, const my::scalar* y);

	void
	scale( const my::real alpha, my::scalar* x);
	my::real
	nrm2( const my::scalar* x);

	void
	normalize( my::scalar* x );

	void
	axpy( 	const my::real alpha,
			const my::scalar* x,
			my::scalar* y);

	void
	swap( my::scalar* x , my::scalar* y);

	void
	copy( my::scalar* x , my::scalar* y);

	private:
		my::integer shift;
		my::integer dim;
		my::integer eff_dim;

	};

#endif /* KPM_KPM_LINALG_HPP_ */
