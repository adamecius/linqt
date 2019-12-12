
#include "chebyshev_moments.hpp"

chebyshev::Moments2D::Moments2D( std::string momfilename )
{
	//Check if the input_momfile have the right extension 
	std::size_t ext_pos = string( momfilename ).find(".chebmom2D"); 
	if( ext_pos == string::npos )
	{ std::cerr<<"The first argument does not seem to be a valid .chebmom2D file"<<std::endl; assert(false);}

	//if it does, use it to get the extension
	system_label = momfilename.substr(0,ext_pos);


	//and then try to open the file
	std::ifstream momfile( momfilename.c_str() );
	assert( momfile.is_open() );

	//if succesful, read the header
	momfile>>this->system_size>>this->band_width>>this->band_center; //in the file what you have is the bandwidth
	momfile>>this->numMoms[0]>>this->numMoms[1];
		
	//create the moment array and read the data
	mu = vector<complex<double> >(numMoms[1]*numMoms[0], 0.0);	
	double rmu,imu;
	for( int m0 = 0 ; m0 < numMoms[0] ; m0++)
	for( int m1 = 0 ; m1 < numMoms[1] ; m1++)
	{ 
		momfile>>rmu>>imu;
		this->operator()(m0,m1) = std::complex<double>(rmu,imu);
	}
	momfile.close();
};
