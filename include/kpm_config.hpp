#ifndef KPM_KPM_CONFIG_HPP
#define KPM_KPM_CONFIG_HPP

#include "types_definitions.hpp"
#include "parser.hpp" //GetValue	

namespace kpm
{
	
	struct config
	{
		std::string label, header, suffix;
		size_t maxMom,maxParMom;
		qt::real energyShift,energyScale, cutOff;

		config(): 
		label(""), header(""), suffix("")
		,maxMom(0),maxParMom(0)
		,energyShift(0.),energyScale(0.), cutOff(0.)
		{};
		
		inline 
		qt::real RescaleToChebDom(const qt::real x )
		{
			return  2.0*cutOff*( x	- energyShift )/energyScale ;
		}

		inline 
		qt::real RescaleToEnergy( qt::real x )
		{
			return  0.5*x*energyScale/cutOff + energyShift  ;
		}
	};


};
	

#endif
