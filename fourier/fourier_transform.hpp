/*
 * fourier_transform.hpp
 *
 *  Created on: Aug 30, 2016
 *      Author: jgarcia
 */

#ifndef FOURIER_TRANSFORM_HPP_
#define FOURIER_TRANSFORM_HPP_

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>

#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>

#include "types_definitions.hpp"
#include "mkl_types.h"
#include "mkl.h"
#include "mkl_service.h"
#include "mkl_dfti.h"

namespace mkl
{

void
direct_FFT(	DFTI_DESCRIPTOR_HANDLE &dft_handle,const my::integer _dim0,
			const my::integer _dim1,
			const my::integer _dim2,
			std::complex<float>* _A);

void
inverse_FFT(DFTI_DESCRIPTOR_HANDLE &dft_handle,const my::integer _dim0,
			const my::integer _dim1,
			const my::integer _dim2,
			std::complex<float>* _A);

};

#endif /* FOURIER_TRANSFORM_HPP_ */
