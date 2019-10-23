
#ifndef LATTICE_GEOMETRY
#define LATTICE_GEOMETRY 

#include <cmath>
#include <complex>

typedef std::complex<double>  Complex;
typedef double  Real;

const static Complex I(0.,1.);
const static int NORB= 2;
const static int SPIN= 2;

// One of the usual lattice vectors	
//const static double A0[]  	 = { 1.0*sqrt(3.), 0.0 ,0 };			// FIRST LATTICE VECTOR
//const static double A1[]  	 = { 0.5*sqrt(3.), 1.5 ,0 };			// SECOND LATTICE VECTOR
//const static double A2[]  	 = { 0.0		 , 0.0 ,1 };			// SECOND LATTICE VECTOR
//const static double Delta[]  = {-0.5*sqrt(3.),-0.5 ,0 };			// THIRD LATTICE VECTOR


// One of the usual lattice vectors	

const static Real 
A[3][3] = { { sqrt(3.) , 0.0 ,0 }, { 0.5*sqrt(3.) ,1.5 ,0 } , { 0.0 , 0.0 ,0 } } ;

const static Real 
Delta[]  = { -0.5*sqrt(3.),-0.5 ,0 };			// THIRD LATTICE VECTOR


const static Complex
SIGMA[3][2][2]= 
{
  { { 0 , 1 },{ 1 , 0 } }	,
  { { 0 ,-I },{ I , 0 } }	,
  { { 1 , 0 },{ 0 ,-1 } }	
};


Complex CrossProductDotZ(const Real x[3], 
						 const int s0, const int s1  )
{
		Real norm= sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ) ; 
		
		return SIGMA[0][s0][s1]*x[1] - SIGMA[1][s0][s1]*x[0];
}




#endif
