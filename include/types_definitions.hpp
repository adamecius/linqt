/*
 * type_config.hpp
 *
 *  Created on: Jun 9, 2016
 *      Author: jgarcia
 */

#ifndef TYPES_DEFINITIONS_HPP_
#define TYPES_DEFINITIONS_HPP_

#ifndef NEIGEN
#include "Eigen/Sparse"
#endif

#ifdef __NVCC__

#include "thrust/complex.h"
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"

typedef float 				  external_real;
typedef thrust::complex<external_real> external_complex;
typedef thrust::complex<external_real > external_scalar;
typedef thrust::host_vector<external_scalar > external_vector;
typedef thrust::device_vector<external_scalar > external_dvector;

#else

#include <complex>
#include <vector>

typedef float 			   external_real;
typedef std::complex< external_real > external_complex;
typedef std::complex< external_real > external_scalar;
typedef std::vector< external_scalar >  external_vector;
typedef std::vector< external_scalar >  external_dvector;

#endif


namespace my
{

typedef external_real 		real;
typedef unsigned long 		size_t;
typedef int					integer;
typedef external_complex	complex;
typedef external_scalar		scalar;
typedef external_vector		vector;
typedef external_dvector	dvector;
#ifndef NEIGEN
typedef Eigen::SparseMatrix< scalar,Eigen::RowMajor, my::integer > SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<scalar> spEntry;
#endif
};

namespace NumCal
{
typedef external_real 		real;
typedef unsigned long 		size_t;
typedef int					integer;
typedef external_complex	complex;
typedef external_scalar		scalar;
typedef external_vector		vector;
typedef external_dvector	dvector;
#ifndef NEIGEN
typedef Eigen::SparseMatrix< scalar,Eigen::RowMajor, my::integer > SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<scalar> spEntry;
#endif

}
const static NumCal::complex I(0.0,1.0);
#endif /* TYPE_CONFIG_HPP_ */
