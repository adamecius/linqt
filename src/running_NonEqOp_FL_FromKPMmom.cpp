
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
const double kb    = 0.086173324; 	// boltlzmann constant meV/K
const double hbar  = 0.6582119624; 	// reduced  planck constant eV.fs

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
	S_OPR = cfg.lookup("Operator");


	int maxNumMom,maxNumTMom,numRV;
	cfg.lookupValue("NumberOfMoments", maxNumMom );
	cfg.lookupValue("NumberOfTMoments", maxNumTMom );
	cfg.lookupValue("NumberOfRandVec", numRV );

	//Use the data read from input file to determine the prefix of all files
	std::string prefix  ="NonEqOp"+S_OPR+LABEL+"KPM_M"+to_string(maxNumMom)+"RV"+to_string(numRV);
	std::string momfilename =prefix+".mom2D" ;


	//Open file, read first the number of moments by reading the
	//index of the last element
	std::cout<<std::endl<<"Reading moments from "<<momfilename<<std::endl;
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
	std::cout<<"Finished moments from "<<momfilename<<std::endl;

	if( numMoms0 != maxNumMom || numMoms1 != maxNumTMom )
	{
		std::cerr	<<std::endl<<"WARNING!"<<std::endl
					<<" The maximum number of moments obtained from the momfile: "
					<<numMoms0<<","<<numMoms1<<std::endl
					<<"do not correspond with the maximum number of moments from the config file: "<<maxNumMom<<","<<maxNumTMom<<std::endl;
		return -1;
	}
	
	//Read the effective number of moments
	int momCutOff = numMoms0,tmomCutOff = numMoms1;
	cfg.lookupValue("MomentsCutOff" ,momCutOff);
	cfg.lookupValue("TMomentCutOff" ,tmomCutOff);
	if (momCutOff <= numMoms0 && momCutOff>0 && tmomCutOff <= numMoms1 && tmomCutOff>0 )
	{
		numMoms0= momCutOff; 	numMoms1= tmomCutOff;
	}	
	else
		std::cerr<<std::endl<<"Warning."<<std::endl
				 <<"Your MomentsCutOff="<<momCutOff<<std::endl
				 <<"or your TMomentsCutOff="<<tmomCutOff<<" is invalid."<<std::endl
				 <<"Using the maximum values= "<<numMoms0<<"x"<<numMoms1<<" as cutoff"<<std::endl;

