
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


struct chebMom
{
	chebMom(const int m0,const int m1):numMom1(m1),numMom0(m0){};
	 
	complex& operator()(const int m0,const int m1)
	{
		return pmu[ m0*numMom1 + m1 ];
	}
	
	int numMom1,numMom0;
	complex* pmu;
};

int main(int argc, char *argv[])
{	
	libconfig::Config cfg;
	// Read the file. If there is an error, report it and exit.
	cfg.readFile(argv[1]);
	
	const
	std::string 
	LABEL = cfg.lookup("SystemName"),
	S_OPR = cfg.lookup("Operator"),
	S_OPL = "VX";

	std::cout<<"Computing the chebyshev moments needed for electrical response of the operator "<<S_OPR<<std::endl;
	std::cout<<"in the system labeled as "<<LABEL<<std::endl;
	
	double bandWidth, bandCenter;
	cfg.lookupValue("BandWidth",bandWidth);
	cfg.lookupValue("BandCenter",bandCenter);
	std::cout<<"The spectrum is considered as extending from "<<0.5*(bandCenter-bandWidth)<<" to "<< 0.5*(bandCenter+bandWidth)<<std::endl;


	int numMoms, numTMoms;
	cfg.lookupValue("NumberOfMoments",numMoms);
	cfg.lookupValue("NumberOfTMoments",numTMoms); if(numTMoms==0) numTMoms=numMoms; 
	std::cout<<"The total number of moments will be "<<numMoms<<"x"<<numTMoms<<std::endl;

	
	int numRV=1, seed=1521;
	cfg.lookupValue("NumberOfRandVec",numRV);
	std::cout<<"Using a total number of random vectors of "<<numRV<<" with seed "<<seed<<std::endl;
	
	std::stringstream ss;
	std::string sNumMom, snumRV;
	ss << numMoms; sNumMom = ss.str(); ss.str(std::string());
	ss << numRV; snumRV = ss.str();

	
	//READ NORMALIZED HAMILTONIAN
	CSRMatrix NHAM,OPL,OPR;
	NHAM.ReadCSRMatrix( "operators/"+LABEL+".HAM.CSR"); //Normalized Hamiltonian
	NHAM.Rescale( 2.0/bandWidth  );
	NHAM.Optimize(numMoms*numMoms*numRV); //OPTIMIZE MATRIX

	//READ OPERATOR LEFT
	OPL.ReadCSRMatrix( "operators/"+LABEL+"."+S_OPL+".CSR");
	OPL.Optimize(numMoms*numRV);

	OPR.ReadCSRMatrix( "operators/"+LABEL+"."+S_OPR+".CSR");
	OPR.Optimize(numRV);
	


	//INITIALIZE KPM PARAMETERS
	const int numMoms0=numMoms,numMoms1=numTMoms; 
	chebMom mu(numMoms0,numMoms1); 
	const int dim = NHAM.Dim();
	complex* data = new complex[ 5*dim + numMoms0*numMoms1 ];
	complex* Psi = &data[0*dim];
	complex* Jm0 = &data[1*dim];
	complex* Jm1 = &data[2*dim];
	complex* Jnm0= &data[3*dim];
	complex* Jnm1= &data[4*dim];
		  mu.pmu = &data[5*dim];
	complex* JT;

	//OUTPUT VARIABLE
	srand(seed);

	//INITIALIZE ITERATION
	for( int r = 0 ; r < numRV ; r++) 
	{
		//DRAW A RANDOM VECTOR
		CreateRandomVector( dim, Jm0 );
		OPR.Multiply(1.0, Jm0, 0.0 , Psi ); // O|Psi>

		//CHEBYSHEV ITERATION
		NHAM.Multiply(1.0, Jm0, 0.0 , Jm1 ); 
		for( int m = 0 ; m < numMoms0; m ++)
		{
			OPL.Multiply(1.0, Jm0, 0.0 , Jnm0 ); // O|Psi>
			NHAM.Multiply(1.0, Jnm0, 0.0 , Jnm1 );
			for( int n = 0 ; n < numMoms1; n ++)
			{
				mu(m,n) += dot(dim ,Psi, Jnm0 );
				JT=Jnm0; Jnm0=Jnm1; Jnm1=JT;
				NHAM.Multiply(2.0, Jnm0, -1.0 , Jnm1 );   // J1 = 2*H*J0
			}
			JT=Jm0; Jm0=Jm1; Jm1=JT;
			NHAM.Multiply(2.0, Jm0, -1.0 , Jm1 );   // J1 = 2*H*J0
		}
	}


	std::ofstream outputfile( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom+"RV"+snumRV+".mom2D").c_str() );
	for( int m = numMoms0-1 ; m >= 0 ; m--)
	for( int n = numMoms1-1 ; n >= 0 ; n--)
	{
		complex mu0 = mu(m,n)/numRV; 
		outputfile<<m<<" "<<n<<" "<<mu0.real()<<" "<<mu0.imag()<<std::endl;
	}
	outputfile<<bandWidth<<" "<<dim<<std::endl;
	outputfile.close();

		
	delete [] data;
	return 0;
}


void LoadKPMParams(std::string config_file)
{
	
}
