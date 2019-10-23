/*
 * fourier_transform.cpp
 *
 *  Created on: Aug 30, 2016
 *      Author: jgarcia
 */


#include "fourier_transform.hpp"


void
mkl::direct_FFT( DFTI_DESCRIPTOR_HANDLE &dft_handle,
				 const my::integer _dim0,
					const my::integer _dim1,
					const my::integer _dim2,
					std::complex<float>* _A)
{
	MKL_LONG l[3]={ _dim0,  _dim1, _dim2 };		//used for the FFT
	MKL_LONG dist=l[0]*l[1]*l[2];
	MKL_LONG status;
	status=DftiCreateDescriptor( &dft_handle, DFTI_SINGLE	, DFTI_COMPLEX , 3, l);
	if( 0 !=  status )
	{
		std::cerr<<"The method direct_FFT couldn't create the descriptor"<<std::endl;
		std::exit(-1);
	}
	//status=DftiSetValue(dft_handle, DFTI_INPUT_DISTANCE, dist);
	//if( 0 !=  status )
	///{
//		std::cerr<<"The method direct_FFT couldn't set the variable DFTI_INPUT_DISTANCE to "<<dist<<std::endl;
//		std::exit(-1);
//	}

	status=DftiCommitDescriptor( dft_handle);
	if( 0 !=  status )
	{
		std::cerr<<"The method direct_FFT couldn't commit the descriptor"<<std::endl;
		std::exit(-1);
	}

	status=DftiComputeForward  ( dft_handle, _A );
	if( 0 !=  status )
	{
		std::cerr<<"The method direct_FFT couldn't perform the direct fourier transform"<<std::endl;
		std::exit(-1);
	}

	//	DftiFreeDescriptor (&dft_handle);

};




void
mkl::inverse_FFT(	 DFTI_DESCRIPTOR_HANDLE &dft_handle,
		const my::integer _dim0,
					const my::integer _dim1,
					const my::integer _dim2,
					std::complex<float>* _A)
{
	MKL_LONG l[3]={ _dim0,  _dim1, _dim2 };		//used for the FFT
	MKL_LONG dist=l[0]*l[1]*l[2];
	MKL_LONG status;
	status=DftiCreateDescriptor( &dft_handle, DFTI_SINGLE	, DFTI_COMPLEX , 3, l);
	if( 0 !=  status )
	{
		std::cerr<<"The method inverse_FFT couldn't create the descriptor"<<std::endl;
		std::exit(-1);
	}
	//status=DftiSetValue(dft_handle, DFTI_INPUT_DISTANCE, dist);
	//if( 0 !=  status )
//	{
//		std::cerr<<"The method inverse_FFT couldn't set the variable DFTI_INPUT_DISTANCE to "<<dist<<std::endl;
//		std::exit(-1);
//	}

	status=DftiCommitDescriptor( dft_handle);
	if( 0 !=  status )
	{
		std::cerr<<"The method inverse_FFT couldn't commit the descriptor"<<std::endl;
		std::exit(-1);
	}

	status=	DftiComputeBackward( dft_handle, _A );
	if( 0 !=  status )
	{
		std::cerr<<"The method inverse_FFT couldn't perform the direct fourier transform"<<std::endl;
		std::exit(-1);
	}
//	DftiFreeDescriptor (&dft_handle);
};

