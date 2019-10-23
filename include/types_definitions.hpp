#ifndef TYPES_DEFINITIONS
#define TYPES_DEFINITIONS

#include <vector>
#include <complex>
#include <limits>


/*!
 *  \addtogroup qt
 *  @{
 */

//! This is the proyect namespace. 
namespace qt
{
	/**
	 * @brief	Use to handle integer scalar values.
	 */
	typedef int 	integer;

	/**
	 * @brief	Use to handle memory localization and sizes.
	 */
	typedef size_t	index;

	/**
	 * @brief	Use to handle real scalar values.
	 */
	typedef double real;

	/**
	 * @brief	Use to handle complex scalar values.
	 */
	typedef std::complex<real> complex;

	/**
	 * @brief Use to handle real vectors, compatibles to C++ Standards.
	 */
	typedef std::vector<real> RealVector ;

	/**
	 * @brief Use to handle complex vectors, compatibles to C++ Standards.
	 */
	typedef std::vector<complex> ComplexVector ;

	/**
	 * @brief Use to handle integer vectors, compatibles to C++ Standards.
	 */
	typedef std::vector<integer> IntegerVector ;

		/**
		 * @brief	Constant used to identify the three-dimensional structure of the program 
		 */
	const integer SPATIAL_DIM=3;


}

/*! @} End of Doxygen Groups*/
#endif
