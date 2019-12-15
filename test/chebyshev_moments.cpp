
// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <cstdlib>		/* for std::cout mostly */
#include <cassert>		/* for std::cout mostly */
#include "chebyshev_moments.hpp"

int main()
{
	int
	dim = 100, numMoms0 = 10, numMoms1 = 20; 
	double BandWidth = 1.0, BandCenter = 2.0;
	chebyshev::Moments2D mu0(numMoms0 , numMoms1); 
	mu0.SystemSize( dim ); 
	mu0.BandWidth( BandWidth); 
	mu0.BandCenter( BandCenter ); 
	for( int m0 = 0 ; m0 < numMoms0 ; m0++ )
	for( int m1 = 0 ; m1 < numMoms1 ; m1++ )
		mu0(m0,m1) = std::complex<double>( rand() , rand() );

	//To test, we create a momfile with the right structure
	std::string momfilename = "test.chebmom2D";
	std::ofstream momfile( momfilename.c_str() );
	momfile<<mu0.SystemSize()<< " "<<mu0.BandWidth()<<" "<<mu0.BandCenter()<<std::endl; 
	momfile<<numMoms0<< " "<<numMoms1<<std::endl;
	for( int m0 = 0 ; m0 < numMoms0 ; m0++ )
	for( int m1 = 0 ; m1 < numMoms1 ; m1++ )
		momfile<<mu0(m0,m1).real()<<" "<<mu0(m0,m1).imag()<<std::endl;
	momfile.close();

	//Read this file 
	chebyshev::Moments2D mu( momfilename.c_str() ); 

	//Compared it
	assert ( mu == mu0 );
	assert ( dim == mu.SystemSize() );
	assert ( BandWidth == mu.BandWidth() );
	assert ( BandCenter == mu.BandCenter() );

return 0;
}
