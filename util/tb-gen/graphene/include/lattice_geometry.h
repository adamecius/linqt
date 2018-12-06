
#ifndef LATTICE_GEOMETRY
#define LATTICE_GEOMETRY 

#define XDIR 0
#define YDIR 1

#include <cmath>
#include <complex>

typedef std::complex<double>  Complex;
typedef double  Real;

const static Complex I(0.,1.);
const static int NORB= 2;
const static int SPIN= 2;
const static int DIMK=NORB*SPIN;	
const static double  ACC = 1.44E-10; 	

// One of the usual lattice vectors	

const  Real 
A[3][3] = 
	 { 
		{ ACC*1.5 ,-ACC*0.8660 ,0 },  
		{ ACC*1.5 , ACC*0.8660 ,0 },
		{ 0.0 ,0.0          ,0.0 }
	  };



const  Real 
Delta[]  = { ( A[0][0] + A[1][0] )/3. ,
			 ( A[0][1] + A[1][1] )/3. ,
		0 
	    };			// THIRD LATTICE VECTOR

//const  Real 
//A[3][3] = 
//	 { 
//		{ 1.5 ,+0.5*sqrt(3) ,0 },  
//		{ 1.5 ,-0.5*sqrt(3) ,0 },
//		{ 0.0 ,0.0          ,1.0 }
//	  };

//const  Real 
//Delta[]  = { A[0][0]/3. + A[1][0]/3. ,
//			 A[0][1]/3. + A[1][1]/3. ,0 };			// THIRD LATTICE VECTOR

const Real 
UnitCellArea	  = std::abs( A[0][0]*A[1][1] - A[0][1]*A[1][0] );

const  Real 
B[2][3] = { 
		{ A[1][1]*2.*M_PI/UnitCellArea,-A[1][0]*2.*M_PI/UnitCellArea,0 }, 
	  	{-A[0][1]*2.*M_PI/UnitCellArea, A[0][0]*2.*M_PI/UnitCellArea,0 }
	   };

const Real 
BZONE	  = std::abs( B[0][0]*B[1][1]- B[0][1]*B[1][0]  );


const Complex
SIGMA[3][2][2]= 
{
  { { 0 , 1 },{ 1 , 0 } }	,
  { { 0 ,-I },{ I , 0 } }	,
  { { 1 , 0 },{ 0 ,-1 } }	
};

//defines the spin matrix
const Complex
Sz[]={
		1,0, 0, 0,//AAdw ABdw AAup ABup
		0,1, 0, 0,//BAdw BBdw BAup BBup//The index compresion is:
		0,0,-1, 0,//AAup ABup AAup ABup//( is*norb+ io )*norb*nspin + 
		0,0, 0,-1 //BAup BBup BAup BBup//( js*norb+ jo )
		};	


inline
Complex 
CrossProductDotZ(const Real x[3], const int s0, const int s1  )
{
		Real norm= sqrt( x[0]*x[0] + x[1]*x[1]  ) ; //normalization vector
		return (SIGMA[0][s0][s1]*x[1] - SIGMA[1][s0][s1]*x[0])/norm;
}

inline
Real
Dot(const int  dim, const Real* x,const Real* y )
{
	Real dot=0;
	for(int i=0;i<dim;i++)
		dot = dot+ x[i]*y[i];

	return dot; 
};

inline
Complex
Dot(const int  dim, const Complex* x,const Complex* y )
{
	Complex dot=0;
	for(int i=0;i<dim;i++)
		dot += std::conj(x[i])*y[i];

	return dot; 
}



#endif
