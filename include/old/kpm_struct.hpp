#ifndef SHARED_KPM_STRUCT_HPP
#define SHARED_KPM_STRUCT_HPP

#include "types_definitions.hpp"
#include "parser.hpp" //GetValue	
#include "kpm_kernel.hpp" //GetValue	

namespace kpm
{
	
	struct option
	{
		std::string label, header, suffix;
		size_t maxMom,maxParMom;
		kpm::real energyShift,energyScale, cutOff;

		option(): 
		label(""), header(""), suffix("")
		,maxMom(0),maxParMom(0)
		,energyShift(0.),energyScale(0.), cutOff(0.)
		{};
		
		inline 
		kpm::real RescaleToChebDom(const kpm::real x )
		{
			return  2.0*cutOff*( x	- energyShift )/energyScale ;
		}

		inline 
		kpm::real RescaleToEnergy( kpm::real x )
		{
			return  0.5*x*energyScale/cutOff + energyShift  ;
		}

	};



	bool GetOptionsFromCFG(std::ifstream& configFile, option& opt)
	{
		bool found=true;
		
		std::string sbuffer,header;
		
		GetValue(configFile,"Label" , opt.label, parser::Optional);
		GetValue(configFile,"Suffix", opt.suffix,parser::Optional);
		GetValue(configFile,"MaxMoments" , opt.maxMom ,parser::Optional,parser::fromInit);		
		GetValue(configFile,"MaxParallMoments" , opt.maxParMom ,parser::Optional,parser::fromInit);		
		
		//Go to header KPM/*
		header="KPM"; 
		found*=FindHeader(configFile, header ) != string::npos;
		found*=GetValue(configFile,"EnergyShift", opt.energyShift);
		found*=GetValue(configFile,"EnergyScale", opt.energyScale);
		found*=GetValue(configFile,"CutOff", opt.cutOff);
		configFile.clear();
		configFile.seekg(0, std::ios::beg);
		return found;
	}
	

	bool PrintOptions(const option opt)
	{
		std::cout<<std::endl<<"KPM_OPTIONS"<<std::endl;
		std::cout<<"Label  = "<<opt.label<<std::endl;
		std::cout<<"Suffix = "<<opt.suffix<<std::endl;
		std::cout<<"EnergyShift = "<<opt.energyShift<<std::endl;
		std::cout<<"EnergyScale = "<<opt.energyScale<<std::endl;
		std::cout<<"CutOff = "<<opt.cutOff<<std::endl;
		if( opt.maxMom!=0 )
			std::cout<<"MaxMom = "<<opt.maxMom<<std::endl;
		if( opt.maxParMom!=0 )
			std::cout<<"MaxParallMoments = "<<opt.maxParMom<<std::endl;

	};


	

};
	

#endif
