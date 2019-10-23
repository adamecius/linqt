
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
	
	libconfig::Config cfg;
	// Read the file. If there is an error, report it and exit.
	
	cfg.readFile(argv[1]);
	
	const
	std::string 
	LABEL = cfg.lookup("SystemName"),
	S_OPR = cfg.lookup("Operator");

	std::string sNumMom, snumRV;
	{
		std::stringstream ss;
		int buff;  cfg.lookupValue("NumberOfMoments", buff );
		ss << buff; sNumMom = ss.str(); ss.str(std::string());
		buff;  cfg.lookupValue("NumberOfRandVec", buff );
		ss << buff; snumRV = ss.str(); ss.str(std::string());
	}
	
	std::string pref_name ="NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom+"RV"+snumRV;
	std::string momfilename =pref_name+".mom2D" ;
	std::cout<<"Reading moment file:"<<momfilename<<std::endl;
	std::ifstream momfile( momfilename.c_str() );

	int m0,m1,numMoms0,numMoms1;
	momfile>>m0>>m1;numMoms0=m0+1; numMoms1=m1+1;
	std::cout<<"The maximum number of moments is:"<<numMoms0<<","<<numMoms1<<std::endl;
 	chebMom mu(numMoms0,numMoms1);
	momfile>>mu(m0,m1).real()>>mu(m0,m1).imag();

	for( m0 = numMoms0-2 ; m0 >= 0 ; m0--)
	for( m1 = numMoms1-2 ; m1 >= 0 ; m1--)
		momfile>>m0>>m1>>mu(m0,m1).real()>>mu(m0,m1).imag();

	double HalfWidth; int dim; 
	momfile>>HalfWidth>>dim; HalfWidth=0.5*HalfWidth;
	momfile.close();


	double Emin, Emax, dE;
	cfg.lookupValue("Emin",Emin);
	cfg.lookupValue("Emax",Emax);
	cfg.lookupValue("dE",dE); 

	
	double tphi_min,tphi_max,dtphi;
	cfg.lookupValue("tphi_max",tphi_max); 

	int M0;
        cfg.lookupValue("M0",M0);  
	if ( M0 != -1 )
	{  numMoms0= M0;
           std::stringstream ss;
           ss << M0; sNumMom = ss.str(); ss.str(std::string());
        }
 	pref_name="NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom+"RV"+snumRV;
       
	const double hbar= 0.6582119624; //eV.fs
	const double br = ( hbar/tphi_max );	//real broadening
	const double nb = ( hbar/tphi_max )/HalfWidth; //normalized
	const double  tmax= 10*tphi_max;//for this particular calculation tmax>>tphi_max
	const double ntmax= tmax*HalfWidth/hbar; //normalized time entering int he coefficients
	std::vector<complex> bcoeff(numMoms1);	//one computes as many coefficients as moments in 1 direction
	std::ofstream outputfile;


	outputfile.open( (pref_name+".FL.NEQ").c_str() );
	for( double E=Emin; E < Emax; E+=dE)
	{
		const double nE  = E/HalfWidth;

		double cond =0;
		for( int m0 = 0 ; m0 < numMoms0 ; m0++)
		for( int m1 = 0 ; m1 < numMoms0 ; m1++)//only use the moments up to nunMoms0 the rest if exists are for time evolution
		{
			complex mu0 = dim*mu(m0,m1)/HalfWidth/HalfWidth/M_PI; 
			cond +=  (mu0*ChebWeigthL(m0,nE,nb) ).imag()*
				  ChebWeigthR(m1,nE,nb).imag(); 
		}
		outputfile<<E<<" "<<cond<<" "<<std::endl;
	}
	outputfile.close();

	outputfile.open( (pref_name+"v2.FL.NEQ").c_str() );
	for( double E=Emin; E < Emax; E+=dE)
	{
		const double nE = E/HalfWidth;
		besselCoeff::BesselCArr(ntmax,nE,1.0,nb,numMoms1,  &bcoeff[0] );

		double cond =0;
		for( int m0 = 0 ; m0 < numMoms0 ; m0++)
		for( int mt = 0 ; mt < numMoms1 ; mt++)
		{
			complex mu0 = dim*mu(m0,mt)/HalfWidth/HalfWidth/M_PI;
			cond  += ChebWeigthR(m0,nE,nb).imag()*( mu0* bcoeff[mt]).imag(); 	
		}	
		outputfile<<E<<" "<<cond<<" "<<std::endl;
	}
	outputfile.close();

	outputfile.open( (pref_name+"eta0.FL.NEQ").c_str() );
	for( double E=Emin; E < Emax; E+=dE)
	{
		const double nE = E/HalfWidth;
		besselCoeff::BesselCArr(ntmax,nE,1.0,0.0,numMoms1,  &bcoeff[0] );

		double cond =0;
		for( int m0 = 0 ; m0 < numMoms0 ; m0++)
		for( int mt = 0 ; mt < numMoms1 ; mt++)
		{
			complex mu0 = dim*mu(m0,mt)/HalfWidth/HalfWidth/M_PI;
			cond  += (mu0*bcoeff[mt]).imag()*ChebWeigthR(m0,nE,nb).imag(); 	
		}	
		outputfile<<E<<" "<<cond<<" "<<std::endl;
	}
	outputfile.close();
	
	return 0;

}
