#include <string>
#include <iostream>
#include <fstream>

#include "types_definitions.hpp"
#include "lattice.hpp"


int main()
{
	const std::string config_filename("example_regular.cfg");

	Lattice honeycomb( config_filename );

	std::cout<<"The label \""<<honeycomb.Label()<<"\""<<std::endl
			 <<"was read from the config file "<<config_filename<<std::endl
			 <<"this data is important, please check if is ok."<<std::endl<<std::endl;

	std::cout<<"The lattice parameters, "
			 <<"read from the config file are:"<<std::endl
			 <<"n0= "  <<honeycomb. SiteNumber(0)<<std::endl
			 <<"n1= "  <<honeycomb. SiteNumber(1)<<std::endl
			 <<"n2= "  <<honeycomb. SiteNumber(2)<<std::endl
			 <<"norb= "<<honeycomb.OrbitalNumber()<<std::endl
			 <<"spin= "<<honeycomb.TotalSpin()<<std::endl
			 <<"total number of orbitals = "<<honeycomb.SiteNumber()<<std::endl;

	std::cout<<"The hoppings are going to be readed from the .hop file"<<std::endl;

	honeycomb.ReadHoppingList();

	honeycomb.HopListToRegMat();

	return 0;
}
