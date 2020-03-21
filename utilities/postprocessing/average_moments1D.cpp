
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <cstdlib>

int main (int argc, char** argv)
{
	if ( argc == 3 ) 
	{
		std::cout<<"This utility must be used with the following arguments NumMOMx <files to average>"<<std::endl;
		return 0;
	}
	int
	NumMom0=std::atoi(argv[1]);
	std::vector< std::complex<double> > mu( NumMom0 ,0 ) ;
	
	//Get the control header, which will be compared to the header in
	//all the moments file to see if we are averaging correctly
	std::string controlHeader;
	std::string output_file( argv[3] );
	std::ifstream moment_file( output_file.c_str() );
	std::getline (moment_file,controlHeader);
	moment_file.close();
	output_file+=".AVERAGE.OUT";
	double avg_num=0;
	for( int i=3; i<argc; i++)
	{
		std::string  header;
		std::string moment_filename( argv[i] );
		//Check if the extension is .mom2D 
		if ( moment_filename.substr(moment_filename.find_last_of(".") + 1) == "mom1D" )
		{ 
			avg_num+=1.;
			std::ifstream moment_file( moment_filename.c_str() );
			std::getline (moment_file,header);
			std::cout<<"Reading file "<<moment_filename<<std::endl;
			if( controlHeader.compare( header) != 0 )
			{
				std::cerr	<<"------------------------------"<<std::endl
							<<"The header: "<<header<<std::endl
							<<"present in file "<<moment_filename<<std::endl
							<<"does not match the control header: "<<controlHeader<<std::endl;
				std::cerr	<<"------------------------------"<<std::endl;
				return 0;
			}
			while( !moment_file.eof()  )
			{
				int m0;
				double muRe,muIm;
				moment_file>>m0>>muRe>>muIm;
				if( m0 < NumMom0  )
					mu[m0 ]+= std::complex<double>(muRe,muIm);
			}
			moment_file.close();
		}else
			std::cout<<"Ignoring file: " <<moment_filename<<" because improper extension"<<std::endl;
	}
	
	std::cout<<"We have averaged "<<avg_num<<" files"<<std::endl;

	std::cout<<"Saving the average file in: "<<output_file<<std::endl;
	std::ofstream moment_output( output_file.c_str() );
	moment_output<<controlHeader<<std::endl;
	for(int m0=0;m0<NumMom0;m0++)
		moment_output<<m0<<" "<<mu[m0 ].real()/avg_num<<" "<<mu[m0 ].imag()/avg_num<<std::endl;
	return 0;

};
