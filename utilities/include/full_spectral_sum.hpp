

#ifndef THRUST_CONDUCTIVITY_SUM_H
#define THRUST_CONDUCTIVITY_SUM_H

#include <thrust/functional.h>
#include <thrust/inner_product.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/complex.h>
#include <fstream>
#include "types_definitions.hpp"

typedef  thrust::complex<my::real> thrust_complex;

namespace utility{

namespace kernel
	{
		
		void 
		JacksonFunction(const my::real m , const my::real M, my::scalar& mu);
		void
		LorentzFunction(const my::real, const my::real lambda, 
						const my::real m , 
						const my::real M, my::scalar& mu);
	}
	namespace sum{

		
	struct SpectralConductivity_binary : public thrust::binary_function<thrust_complex, my::integer, thrust_complex >
	{
		const my::integer M0,M1;
		const thrust_complex I;

		my::real En;

		SpectralConductivity_binary(const my::integer _M0,const my::integer _M1):
			M0(_M0) ,M1(_M1) ,  I(0,1)
			{
			}

		__host__ __device__ 
		my::scalar
		operator()( const thrust_complex &mu, const my::integer &mn) const
		{
			const my::integer m=mn/M0;
			const my::integer n=mn%M1;
			const my::real theta_m= m*acos(En);
			const my::real theta_n= n*acos(En);

			return mu*(
						my::real(cos( theta_m ))*thrust::exp(+I*theta_n)*( En-I*my::real(n*sqrt(1-En*En)) ) +
						my::real(cos( theta_n ))*thrust::exp(-I*theta_m)*( En+I*my::real(m*sqrt(1-En*En)))
					  );
		}
	}; // end plus

	void 
	SpectralConductivity(const std::string moment_filename ,const my::integer M0, const my::integer M1,
						const my::integer NE,
						const std::string output_filename);

	};
};

#endif
