#include "chebyshev_moments.hpp"
#include "mkl.h"

chebyshev::Moments2D::Moments2D( std::string momfilename )
{
	//Check if the input_momfile have the right extension 
	std::size_t ext_pos = string( momfilename ).find(".chebmom2D"); 
	if( ext_pos == string::npos )
	{ std::cerr<<"The first argument does not seem to be a valid .chebmom2D file"<<std::endl; assert(false);}

	//if it does, use it to get the extension
	this->SystemLabel( momfilename.substr(0,ext_pos) ); 

	//and then try to open the file
	std::ifstream momfile( momfilename.c_str() );
	assert( momfile.is_open() );

	//if succesful, read the header
	int ibuff; double dbuff;
	momfile>>ibuff; this->SystemSize(ibuff);
	momfile>>dbuff; this->BandWidth(dbuff);
	momfile>>dbuff; this->BandCenter(dbuff);

	//create the moment array and read the data
	
	momfile>>this->numMoms[0]>>this->numMoms[1];

	this->MomentVector( Moments::vector_t(numMoms[1]*numMoms[0], 0.0) );
	double rmu,imu;
	for( int m0 = 0 ; m0 < numMoms[0] ; m0++)
	for( int m1 = 0 ; m1 < numMoms[1] ; m1++)
	{ 
		momfile>>rmu>>imu;
		this->operator()(m0,m1) = Moments::value_t(rmu,imu);
	}
	momfile.close();
};



void chebyshev::Moments2D::saveIn(std::string filename)
{

    typedef std::numeric_limits<double> dbl;
	ofstream outputfile(filename.c_str());
    outputfile.precision(dbl::digits10);
    outputfile << this->SystemSize() << " " << this->BandWidth() << "  " << this->BandCenter() << std::endl;
    //Print the number of moments for all directions in a line
    for ( auto x : numMoms )
        outputfile << x << " ";
    outputfile << std::endl;

    for ( auto mom : this->MomentVector() )
        outputfile << mom.real() << " " << mom.imag() << std::endl;
    outputfile.close();
};


void chebyshev::Moments2D::Print()
{
	std::cout<<"\n\nCHEBYSHEV MOMENTS INFO"<<std::endl;
	std::cout<<"\tSYSTEM:\t\t\t"<<this->SystemLabel()<<std::endl;
	if( this-> SystemSize() > 0 )
		std::cout<<"\tSIZE:\t\t\t"<<this-> SystemSize()<<std::endl;

	std::cout<<"\tMOMENTS SIZE:\t\t"<<"("<<this->HighestMomentNumber(0)<<" x " <<this->HighestMomentNumber(1)<<")"<<std::endl;
	std::cout<<"\tSCALE FACTOR:\t\t"<<this->ScaleFactor()<<std::endl;
	std::cout<<"\tSHIFT FACTOR:\t\t"<<this->ShiftFactor()<<std::endl;
	std::cout<<"\tENERGY SPECTRUM:\t("
			 <<-this->HalfWidth()+this->BandCenter()<<" , "
			 << this->HalfWidth()+this->BandCenter()<<")"<<std::endl<<std::endl;

};
