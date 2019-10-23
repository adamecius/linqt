
#ifndef LATTICE_GEOMETRY
#define LATTICE_GEOMETRY 

#include "types_definitions.hpp"
const static int NORB= 2;
const static int SPIN= 2;

// One of the usual lattice vectors	
//const static double A0[]  	 = { 1.0*sqrt(3.), 0.0 ,0 };			// FIRST LATTICE VECTOR
//const static double A1[]  	 = { 0.5*sqrt(3.), 1.5 ,0 };			// SECOND LATTICE VECTOR
//const static double A2[]  	 = { 0.0		 , 0.0 ,1 };			// SECOND LATTICE VECTOR
//const static double Delta[]  = {-0.5*sqrt(3.),-0.5 ,0 };			// THIRD LATTICE VECTOR


// One of the usual lattice vectors	
const static my::real A[3][3] = { { sqrt(3.) , 0.0 ,0 }, { 0.5*sqrt(3.) ,1.5 ,0 } , { 0.0 , 0.0 ,0 } } ;

//const static my::real Delta[]  = { ( A[0][0]-2* A[1][0])/3. , (  A[0][1]  -2*A[1][1])/3. , 0 };			// THIRD LATTICE VECTOR
const static my::real Delta[]  = { (A[0][0]+A[1][0])/3. , (A[0][1]+A[1][1])/3. , 0 };			// THIRD LATTICE VECTOR
//const static my::real Delta[]  = { (-2*A[0][0]+ A[1][0])/3. , (-2*A[0][1]+A[1][1])/3. , 0 };			// THIRD LATTICE VECTOR


const static double AREA	  = A[0][0]*A[1][1] - A[0][1]*A[1][0] ;





//comparison with all code
//const static my::real A[3][3] =  { { sqrt(3.) , 0.0 ,0 }, { 0.5*sqrt(3.) ,1.5 ,0 } , { 0.0 , 0.0 ,0 } } ;
//const static my::real Delta[]  = { - sqrt(3.)*0.5  ,-0.5 , 0 };			// THIRD LATTICE VECTOR


//const static double B0[]	  ={ A1[1]*2.*M_PI/Area, -A1[0]*2.*M_PI/Area,0};
//const static double B1[]	  ={-A0[1]*2.*M_PI/Area,  A0[0]*2.*M_PI/Area,0};
//const static double B2[]	  ={ 0.0				 , 0. 				,2.*M_PI/A2[2] };



#endif
