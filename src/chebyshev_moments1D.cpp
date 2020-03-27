#include "chebyshev_moments.hpp"

void chebyshev::Moments1D::Print()
{
	std::cout<<"\n\nCHEBYSHEV 1D MOMENTS INFO"<<std::endl;
	std::cout<<"\tSYSTEM:\t\t\t"<<this->SystemLabel()<<std::endl;
	if( this-> SystemSize() > 0 )
		std::cout<<"\tSIZE:\t\t\t"<<this-> SystemSize()<<std::endl;

	std::cout<<"\tMOMENTS SIZE:\t\t"<<this->HighestMomentNumber()<<std::endl;
	std::cout<<"\tSCALE FACTOR:\t\t"<<this->ScaleFactor()<<std::endl;
	std::cout<<"\tSHIFT FACTOR:\t\t"<<this->ShiftFactor()<<std::endl;
	std::cout<<"\tENERGY SPECTRUM:\t("
			 <<-this->HalfWidth()+this->BandCenter()<<" , "
			 << this->HalfWidth()+this->BandCenter()<<")"<<std::endl<<std::endl;

};

	
