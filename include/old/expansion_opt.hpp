#ifndef SHARED_EXPANSION_OPTION_HPP
#define SHARED_EXPANSION_OPTION_HPP

#include "types_definitions.hpp"
#include "parser.hpp" //GetValue	
#include "kpm_kernel.hpp" //GetValue	


namespace expansion
{
struct option
{
	std::string expType, kerType;
	kpm::real broadening;
	int numMom;
	
	void PrintOptions()
	{
			std::cout<<std::endl;
			std::cout<<"EXPANSION_TYPE: "<<expType<<std::endl;
			std::cout<<"Kernel type: "<<kerType<<std::endl;
			std::cout<<"Broadening: "<<broadening<<" meV" <<std::endl;
			std::cout<<"Number of moments: "<<numMom<<" *"<< std::endl;
			std::cout<<"*If submitted, then it preceed Broadening, if zero, then it is defined automatically through broadening"<<std::endl;	
	};
};
}
#endif
