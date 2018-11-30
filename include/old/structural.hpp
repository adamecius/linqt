

#ifndef STRUCTURAL_HPP
#define STRUCTURAL_HPP



//C libraries
#include <cstdlib>
#include <cassert>
//C++99 libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
//Custom Libraries
#include "types_definitions.hpp"
#include "parser.hpp"


struct Structural
{
	Structural()
	{
		spatialDim=3;
		dim[0]=1; dim[1]=1; dim[2]=1;
		norb = 1;
		nspin = 1;
		
	}
	
	size_t spatialDim;
	size_t dim[3];
	kpm::real lat[3][3];
	size_t norb;
	size_t nspin; 
		
		
	kpm::real volume()
	{
		return 	dim[2]*dim[1]*dim[0]*(
				lat[2][0]*( lat[0][1]*lat[1][2] - lat[0][2]*lat[1][1] )+
				lat[2][1]*( lat[0][2]*lat[1][0] - lat[0][0]*lat[1][2] )+
				lat[2][2]*( lat[0][0]*lat[1][1] - lat[0][1]*lat[1][0] ) 
									);
	};
	
	void savetxt(std::string filename)
	{
		std::ofstream os(filename.c_str(),std::ofstream::binary );	
		for(size_t i=0; i< spatialDim; i++)
		{
			os<<dim[i]<<" ";
			for(size_t j=0; j< spatialDim; j++)
				os<<lat[i][j]<<" ";
		}
		os<<norb<<" ";
		os<<nspin<<" ";
		os.close();
	};

	void loadtxt(std::string filename)
	{
		if ( !fileExists(filename) )
		{
			std::cout<<"Structural file: "<<filename<<" not found. ABORTING "<<std::endl;
			std::exit(-1);
		}
		std::ifstream is(filename.c_str(),std::ifstream::binary );	
		for(size_t i=0; i< spatialDim; i++)
		{
			is>>dim[i];
			for(size_t j=0; j< spatialDim; j++)
				is>>lat[i][j];
		}
		is>>norb;
		is>>nspin;
		is.close();
	}

	void GetFromCFG(std::string filename)
	{
		if ( !fileExists(filename) )
		{
			std::cout<<"the .cfg file: "<< filename<<" not found. ABORTING "<<std::endl;
			std::exit(-1);
		}
		std::ifstream is(filename.c_str(),std::ifstream::binary );	

		
		const std::string header("STRUCTURAL");
		std::string line;
		bool header_found=false;
		while( getline(is, line ) ) 
		{
			if (line.find(header, 0) != std::string::npos)
			{
				for(size_t i=0; i< spatialDim; i++)
				{
					is>>dim[i];
				for(size_t j=0; j< spatialDim; j++)
					is>>lat[i][j];
				}
				is>>norb;
				is>>nspin;
			header_found=true;
			}
		}	
		is.close();

		if( !header_found )
		{
			std::cout << "The header " << header<<" was not found in CFG file " <<filename<<" , therefore aborting"<< std::endl;
			std::exit(-1);
		}
	}

	
	void printParams()
	{
	std::cout<<std::endl<<"The structural parameters are:"<<std::endl
			 <<"Unit cell lattice vectors :"<<std::endl
			 <<"Lat1: ( "<<lat[0][0]<<" , "<<lat[0][1]<<" , "<<lat[0][2]<<" )"<<std::endl
			 <<"Lat2: ( "<<lat[1][0]<<" , "<<lat[1][1]<<" , "<<lat[1][2]<<" )"<<std::endl
			 <<"Lat3: ( "<<lat[2][0]<<" , "<<lat[2][1]<<" , "<<lat[2][2]<<" )"<<std::endl
			 <<"Supercell with dimensions : "
			 <<dim[0]<<" x "<<dim[1]<<" x "<<dim[2]<<std::endl
			 <<"Number of orbital per unit cell is : "
			 <<norb<<std::endl
			 <<"Number os spin per orbital spin : "
			 <<nspin
			 <<std::endl;
	}
	
};	  


#endif
