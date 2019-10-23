
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
#include "kubo_bastin_aux.hpp"

const int TIME_NOT_FOUND = -2 ;
const int INFINITE_TIME = -1 ;
const int INTEGER_NOT_FOUND = 0 ;
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
 	chebMom ChebMom(numMoms0,numMoms1);
	momfile>>ChebMom(m0,m1).real()>>ChebMom(m0,m1).imag();

	for( m0 = numMoms0-2 ; m0 >= 0 ; m0--)	//energy index
	for( m1 = numMoms1-2 ; m1 >= 0 ; m1--)	//time index
		momfile>>m0>>m1>>ChebMom(m0,m1).real()>>ChebMom(m0,m1).imag();


	double HalfWidth; int dim; 
	momfile>>HalfWidth>>dim; HalfWidth=0.5*HalfWidth;
	momfile.close();
	
	double Emin, Emax, dE;
	cfg.lookupValue("Emin",Emin);
	cfg.lookupValue("Emax",Emax);
	cfg.lookupValue("dE",dE); 


	double tphi=TIME_NOT_FOUND;	//Default value. when not found keep the value same
	cfg.lookupValue("tphi",tphi); 


	double tmax=TIME_NOT_FOUND;	//Default value. when not found keep the value same
	cfg.lookupValue("tmax",tmax); 
	
	int M0 = INTEGER_NOT_FOUND ;
	cfg.lookupValue("MomentsCutoff",M0);  
	if ( M0 != INTEGER_NOT_FOUND )
	{  
		numMoms0= M0;
		numMoms1= M0;
		std::stringstream ss;
		ss << M0; sNumMom = ss.str(); ss.str(std::string());
    }

	double temp = 0;
	cfg.lookupValue("temperature",temp); 

	const double time_scaling = 8;
	const double hbar = 0.6582119624; //eV.fs
	const double kb   = 8.6173324E-5; //eV/K
	const double ntemp= temp*kb/HalfWidth; //normalized time entering int he coefficients


	//Set an automatic dephasing if not define
	double br;
	if( tphi == TIME_NOT_FOUND )	//Then use an automatic maximum time
		tphi = time_scaling*hbar/temp/kb;	//
	if( tphi == INFINITE_TIME )
		br = 0.0;
	else
		br = ( hbar/tphi );				//real broadening
	const double nb   = ( hbar/tphi )/HalfWidth; 	//normalized


	//Set an automatic time if not define
	if( tmax == TIME_NOT_FOUND )	//Then use an automatic maximum time
		tmax = time_scaling*tphi;//for this particular calculation tmax>>tphi_max
	const double ntmax = tmax*HalfWidth/hbar; //normalized time entering int he coefficients


	std::vector<complex> bcoeff(numMoms1);	//one computes as many coefficients as moments in 1 direction
	std::ofstream outputfile;

	std::cout<<"The calculation will be performed using the following set of parameters:"<<std::endl
			 <<"Temperature: "<<temp<<" K "<<" or beta: "<<hbar/temp/kb<<" (fs) "<<std::endl
			 <<"dephasing time: "<<tphi<<" fs or broadening "<<1000*br<<" meV "<<std::endl;
	if( tmax == INFINITE_TIME )
		std::cout<<"evolution time: infinity  "<<std::endl;
	else 
		std::cout<<"evolution time: "<<tmax<<" fs "<<std::endl;

	std::cout<<"Chebyshev Moments: ("<<numMoms0<<","<<numMoms1<<")"<<std::endl;

	//Initialize a ser of alpha coefficients
	KuboBastin::AlphaCoeff  alpha_mn(numMoms0, numMoms1);
	//Compute the intermediate GammaCoefficients

	outputfile.open( (pref_name+"KuboGreenwood.COND").c_str() );
	std::cout<<"Computing the NonEquilibrium Kubo Greenwood formula for "<<ntemp<<std::endl;
	{
		alpha_mn.ComputeGammaCoeff_FL(ntmax, nb );
		const int nim = int( (Emax-Emin)/dE );
		const double units = dim/HalfWidth/HalfWidth;
		for(int im= 0 ; im < nim; im++ )
		{
			double mu = Emin + (Emax-Emin)*im/(double)(nim-1);
			const double nmu = mu/HalfWidth;
			alpha_mn.ComputeAlphaCoeff_FL( nmu, ntemp);
			const complex I(0,1);
			double output = 0.0;
			for( int m1 = 0 ; m1 < numMoms1; m1++)
			for( int m0 = 0 ; m0 < numMoms0; m0++)
			{
				double scal = 2.0;
				if( m0 == 0 ) scal*=0.5;
				if( m1 == 0 ) scal*=0.5;
				output += 1.0*(alpha_mn(m0,m1)*ChebMom(m0,m1) ).imag()*scal;//the two is due to the imaganry part
			}
			outputfile<<mu<<" "<<output*units<<std::endl;
		}
	}
	outputfile.close();
	std::cout<<"Finish. Ending Calculation"<<ntemp<<std::endl;

	outputfile.open( (pref_name+"KuboBastin.COND").c_str() );
	std::cout<<"Computing the NonEquilibrium Kubo Bastin formula for "<<ntemp<<std::endl;
	{
		alpha_mn.ComputeGammaCoeff(ntmax, nb );
		const int nim = int( (Emax-Emin)/dE );
		const double units = dim/HalfWidth/HalfWidth;
		for(int im= 0 ; im < nim; im++ )
		{
			double mu = Emin + (Emax-Emin)*im/(double)(nim-1);
			const double nmu = mu/HalfWidth;
			const complex I(0,1);
			double output = 0.0;
			for(double  dir= -1.0 ; dir <= 1.0; dir += 2.0 )
			{ 	//This define the ingration direction form right to lef (-1) or viceversa (1)
				alpha_mn.ComputeAlphaCoeff( nmu, ntemp,dir);
				for( int m1 = 0 ; m1 < numMoms1; m1++)
				for( int m0 = 0 ; m0 < numMoms0; m0++)
				{
					double scal = 2.0;
					if( m0 == 0 ) scal*=0.5;
					if( m1 == 0 ) scal*=0.5;
					output += -1.0*(alpha_mn(m0,m1)*ChebMom(m0,m1) ).imag()*scal; //the two is due to the imaganry part is cancel by the 2 due to average
				}
			}
			outputfile<<mu<<" "<<output*units<<std::endl;
		}
	}
	outputfile.close();

	return 0;

}
