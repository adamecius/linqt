
// Used for OPENMP functions
#include <omp.h>

// Parsing library
#include <libconfig.h++>

// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/

// MKL LIBRARIES
#define MKL_Complex16 std::complex<double>
#include "mathimf.h"
#include "mkl.h"
#include "mkl_spblas.h"
typedef std::complex<double> complex;
typedef MKL_INT integer;


//CUSTOM LIBRARIES
#include "kernel_functions.hpp"
#include "algebra_functions.hpp"
#include "bessel_int_coeff.hpp"
#include "time_handler.hpp"

const int INTEGER_NOT_FOUND = 0 ;

template <typename T> 
std::string to_string( T x ) 
{
	std::stringstream ss;
	ss << x;
	return ss.str();
}


struct chebMom
{
	chebMom():numMom1(0),numMom0(0){};

	chebMom(const int m0,const int m1):numMom1(m1),numMom0(m0),mu( std::vector<complex>(numMom1*numMom0) ){};
	 
	complex& operator()(const int m0,const int m1)
	{
		return mu[ m0*numMom1 + m1 ];
	}
	int numMom1,numMom0;
	std::vector<complex> mu;
};


int main(int argc, char *argv[])
{	
	libconfig::Config cfg;	//Class use to handle files
	cfg.readFile(argv[1]); 	// Read the file. If there is an error, report it and exit.
	
	const
	std::string 
	LABEL = cfg.lookup("SystemName"),
	S_OPR = "VX";

	int maxNumMom, numRV;
	cfg.lookupValue("NumberOfMoments", maxNumMom );
	cfg.lookupValue("NumberOfRandVec", numRV );

	std::cout<<"Computing the conductivity for the operator "<<S_OPR<<std::endl
			 <<"Using a maximum of "<<maxNumMom<<" chebyshev moments"<<std::endl
			 <<"computed using "<<numRV<<" random vectors"<<std::endl;
	

	//Use the data read from input file to determine the prefix of all files
	std::string prefix  ="NonEqOp"+S_OPR+LABEL+"KPM_M"+to_string(maxNumMom)+"RV"+to_string(numRV);
	std::string momfilename =prefix+".mom2D" ;
	std::cout<<"Reading moments from "<<momfilename<<std::endl;


	//Open file, read first the number of moments by reading the
	//index of the last element
	std::ifstream momfile( momfilename.c_str() );
	int m0,m1,numMoms0,numMoms1;
	momfile>>m0>>m1;numMoms0=m0+1; numMoms1=m1+1;
	//Save these moments into a chebMom file
 	chebMom mu(numMoms0,numMoms1);
	momfile>>mu(m0,m1).real()>>mu(m0,m1).imag();
	for( m0 = numMoms0-2 ; m0 >= 0 ; m0--)	//energy index
	for( m1 = numMoms1-2 ; m1 >= 0 ; m1--)	//time index
		momfile>>m0>>m1>>mu(m0,m1).real()>>mu(m0,m1).imag();
	//and finally read the halfwidth and the dimension at the end of the
	//moment file
	double HalfWidth; int dim; 
	momfile>>HalfWidth>>dim; HalfWidth=0.5*HalfWidth;
	momfile.close();

	if( numMoms0 != maxNumMom || numMoms1 != maxNumMom )
	{
		std::cerr	<<"WARNING!"<<std::endl
					<<" The maximum number of moments obtained from the momfile: "
					<<numMoms0<<","<<numMoms1<<std::endl
					<<"do not correspond with the maximum number of moments from the config file: "<<maxNumMom<<std::endl;
		return -1;
	}

	//Read fron config file the energy grid where the conductivity will be computed
	double Emin=-1.0, Emax=1.0, dE=0.1;
	try{
		const libconfig::Setting &enerGrid = cfg.lookup("EnergyGrid");
		std::vector<double> enerGridParam( enerGrid.getLength() );
		for(int i = 0; i < enerGrid.getLength(); ++i)
			enerGridParam[i] = (double)enerGrid [i];
		Emin= enerGridParam[0]; 	
		Emax= enerGridParam[1];
		dE  = enerGridParam[2];
		}
	catch(const libconfig::SettingNotFoundException &nfex)
	{
		std::cerr<<"EnergyGrid field not found, using default values"<<std::endl;
	}
	std::cout	<<"The conductivity will be computed in the uniform grid"<<std::endl
				<<"between "<<Emin<<" and "<<Emax<<" with steps "<<dE<<" in units of eV"<<std::endl;

	//Read the temperature, dephasing time, and simulation time 
	const double kb    = 0.086173324; 	// boltlzmann constant meV/K
	const double hbar  = 0.6582119624; 	// reduced  planck constant eV.fs

	//TEMPERATURE
	double temp=300.0; 
	cfg.lookupValue("Temperature",temp); 

	//TIMES ARE HANDLED USING TIME UNIT TO HANDLE INFINITIES
	time_unit tphi = 10.0*temp*kb/hbar ;
	tphi.ReadTime(cfg,"tphi" );
	time_unit tmax; if( !tphi.Infinite() ) tmax = 10.*tphi ; else tmax = 10.0*temp*kb/hbar;
	tmax.ReadTime(cfg,"tmax" );

	//Scaling  and normalization
	const double tbr   	= temp*kb;		//thermal broadening
	const double ntb   	= tbr/HalfWidth; 	//normalized thermal broadening
	time_unit ntphi 	= tphi*HalfWidth/hbar; //normalized tphi
	time_unit ntmax 	= tmax*HalfWidth/hbar; //normalized tmax

	//Read the effective number of moments
	int momentCutOff = maxNumMom;
	cfg.lookupValue("momentCutOff",momentCutOff);
	if (momentCutOff > maxNumMom || momentCutOff<0 )
	{
		std::cerr<<"Warning."<<std::endl
				 <<"You are requiring more moments"<<momentCutOff<<" than the maximum"<<maxNumMom<<std::endl
				 <<"using the maximum value as cutoff"<<std::endl;
		momentCutOff = maxNumMom ;
    }
	numMoms0= momentCutOff; 	numMoms1= momentCutOff;

	if( tphi.Infinite() &&tmax.Infinite() )
	{
		std::cerr<<" tphi and tmax cannot be both infinite simultaneously"<<std::endl;
		return -1; 
	}
	

	std::cout	<<"The calculation will be performed using: "<<std::endl
				<<" temperature = "<<temp<<" K corresponding to "<<hbar/tbr<<" fs or "<<tbr<<" meV"<<std::endl
				<<" tphi = "<<tphi<<" fs corresponding to "<<1000*hbar/tphi<<" meV"<<std::endl
				<<" tmax = "<<tmax<<" fs corresponding to "<<1000*hbar/tmax<<" meV"<<std::endl
				<<" numberOfMoments = "<<numMoms0<<" x "<<numMoms1<<std::endl;


	std::string
	outputName  ="COND"+LABEL+"KPM_M"+to_string(numMoms0)+"x"+to_string(numMoms1)+"RV"+to_string(numRV)+"tphi"+(std::string)tphi+"tmax"+(std::string)tmax+".dat";
	std::cout<<"Saving the data in "<<outputName<<std::endl;
	
	std::vector<complex> bcoeff(  numMoms1 );
	std::ofstream outputfile( outputName.c_str() );
	for( double E=Emin; E < Emax; E+=dE)
	{
		const double nE = E/HalfWidth;
		besselCoeff::BesselCArr(ntphi,0.0,ntmax,nE,numMoms1,&bcoeff[0] );
		double cond =0;
		for( int m0 = 0 ; m0 < numMoms0 ; m0++)
		for( int mt = 0 ; mt < numMoms1 ; mt++)
		{
			double mu0 = dim*mu(m0,mt).real()/HalfWidth/HalfWidth/M_PI;
			cond  += mu0*ChebWeigthR(m0,nE,ntb).imag()*bcoeff[mt].imag(); 	
		}	
		outputfile<<E<<" "<<cond<<" "<<std::endl;
	}	
	outputfile.close();


	return 0;

}
