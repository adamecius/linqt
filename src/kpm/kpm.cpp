#include "kpm.hpp"


bool
kpm::GetOptionsFromCFG(std::ifstream& configFile, config& opt)
{
	bool found=true;

	std::string sbuffer,header;

	parser::GetValue(configFile,"Label" , opt.label, parser::Optional);
	parser::GetValue(configFile,"Suffix", opt.suffix,parser::Optional);
	parser::GetValue(configFile,"MaxMoments" , opt.maxMom ,parser::Optional,parser::fromInit);

	//Go to header KPM/*
	header="KPM";
	found*=parser::FindHeader(configFile, header ) != std::string::npos;
	found*=parser::GetValue(configFile,"EnergyShift", opt.energyShift);
	found*=parser::GetValue(configFile,"EnergyScale", opt.energyScale);
	found*=parser::GetValue(configFile,"CutOff", opt.cutOff);
	configFile.clear();
	configFile.seekg(0, std::ios::beg);
	return found;
}

void
kpm::PrintOptions(const config opt)
{
	std::cout<<std::endl<<"KPM_OPTIONS"<<std::endl;
	std::cout<<"Label  = "<<opt.label<<std::endl;
	std::cout<<"Suffix = "<<opt.suffix<<std::endl;
	std::cout<<"EnergyShift = "<<opt.energyShift<<std::endl;
	std::cout<<"EnergyScale = "<<opt.energyScale<<" (should be >0 )"<<std::endl;
	std::cout<<"CutOff = "<<opt.cutOff<<" (should be >0 )"<<std::endl;
	std::cout<<"MaxMom = "<<opt.maxMom<<" (should be >0 )"<<std::endl;
};


void kpm::cheb_evolve(	const qt::index n, 
						const sparse::BlockMatrix<qt::complex>& H,
						qt::complex* &jn0,
						qt::complex* &jn1
						)
{
	switch (n)
	{
		case 0:
		{
			H.Multiply( 1. , jn0 , 0 , jn1 );
			break;
		}

		default:
		{
			qt::complex *jnt;
			H.Multiply( 2. , jn1 ,-1. , jn0 );
			jnt= jn1 ;
			jn1= jn0;
			jn0= jnt;
			break;
		}
	}
}

