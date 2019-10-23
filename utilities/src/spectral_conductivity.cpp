
#include "full_spectral_sum.hpp"

int main (int argc, char** argv)
{

	const int M0=atoi(argv[1]);	

	if(M0!=0)
		std::cout<<"You selected the number of moments"<<M0<<std::endl;
	else
	std::cout<<"Using the total number of moments specified in the file"<<std::endl;
	for( int i=2; i<argc; i++)
	{
		std::string inputxx( argv[i] );
		std::string outputxx("conductivity"+inputxx+".dat");
	        utility::sum::SpectralConductivity(inputxx,M0,M0, 55001,outputxx);

	}
 


	return 0;

};
