#ifndef TYPES_DEFINITIONS
#define TYPES_DEFINITIONS

#include <complex>
#include <limits>


namespace qt
{

typedef double real;
typedef std::complex<real> complex;
typedef long int 	integer;
typedef long int	index;
typedef long int 	dimension;

const qt::complex I(0.,1.);
}

#endif
