
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
#include <limits>		//For getting machine precision limit

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



int estimate_batch_size(const int numMoms)
{
	
}



template <typename T> 
std::string to_string( T x ) 
{
	std::stringstream ss;
	ss << x;
	return ss.str();
}
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
	
	double HalfWidth=10, bandCenter=0;
	cfg.lookupValue("BandWidth",HalfWidth);HalfWidth=HalfWidth/2.0;
	cfg.lookupValue("BandCenter",bandCenter);
	std::cout<<"The spectrum is considered as extending from "<<(bandCenter-HalfWidth)<<" to "<< (bandCenter+HalfWidth)<<std::endl;


	int numMoms=0, numTMoms=0;
	cfg.lookupValue("NumberOfMoments",numMoms);
	cfg.lookupValue("NumberOfTMoments",numTMoms); if(numTMoms==0) numTMoms=numMoms; 
	std::cout<<"The total number of moments will be "<<numMoms<<"x"<<numTMoms<<std::endl;
	std::string
	sNumMom0 = to_string(numMoms), 
	sNumMom1 = to_string(numTMoms);


	
	//READ NORMALIZED HAMILTONIAN
	CSRMatrix NHAM,OPL,OPR;
	std::cout<<"Reading operator "+LABEL+".HAM.CSR"<<std::endl;
	NHAM.ReadCSRMatrix( "operators/"+LABEL+".HAM.CSR"); //Normalized Hamiltonian
	NHAM.Rescale( 1.0/HalfWidth  );
	//NHAM.Rescale( 2.0/bandWidth  );
	//READ OPERATOR LEFT
	std::cout<<"Reading operator "+LABEL+"."+S_OPL+".CSR"<<std::endl;
	OPL.ReadCSRMatrix( "operators/"+LABEL+"."+S_OPL+".CSR");
	//READ OPERATOR RIGHT
	std::cout<<"Reading operator "+LABEL+"."+S_OPR+".CSR"<<std::endl;
	OPR.ReadCSRMatrix( "operators/"+LABEL+"."+S_OPR+".CSR");

	int numSV=1, blockSize = NHAM.Dim(), blockShift = 1, seed = 1541581, bvecStride=1 , bvecSize =blockSize ;  
	int numBlock = (NHAM.Dim() + blockSize -1)/ blockSize ;
	cfg.lookupValue("NumberOfStateVec",numSV);
	cfg.lookupValue("NumberOfBlocks",numBlock);
	cfg.lookupValue("StateBlockSize",blockSize);
	cfg.lookupValue("StateBlockShift",blockShift);
	cfg.lookupValue("BlockVectorStride",bvecStride);
	cfg.lookupValue("BlockVectorSize",bvecSize);  
	std::cout<<"Using "<<numSV<<" state vectors in "<< numBlock<<" blocks of size "<<blockSize<<" shifted by "<<blockShift<<" "<<numBlock<<std::endl;
	std::string snumRV=to_string(numSV);

	OPR.Optimize(numSV);
	OPL.Optimize(numMoms*numSV);
	NHAM.Optimize(numMoms*numMoms*numSV); //OPTIMIZE MATRIX

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
	srand (time(NULL));

	//INITIALIZE ITERATION
	for( int r = 0 ; r < numSV ; r++) 
	for( int b = 0 ; b < numBlock ; b++) 
	{
		std::cout<<b <<std::endl;
		const int i0_shift = b*blockSize;
		if( i0_shift < dim )
		{
			//DRAW A RANDOM VECTOR
			for( int i =0; i < dim ; i++ )
				Jm0[i]=0.0;
			std::cout<<"State vector is  running from: "<<i0_shift<<"->"<<bvecSize<<" strided by "<<bvecStride<<" and the system dimension is "<<dim<<std::endl;
			CreateRandomVector( bvecSize,i0_shift, bvecStride	, Jm0 );

			OPR.Multiply(1.0, Jm0, 0.0 , Psi ); // O|Psi>
			//CHEBYSHEV ITERATION
			NHAM.Multiply(1.0, Jm0, 0.0 , Jm1 );
			axpy(dim ,-bandCenter/HalfWidth,Jm0, Jm1 ); 
			for( int m = 0 ; m < numMoms0; m ++)
			{
				OPL.Multiply(1.0, Jm0, 0.0 , Jnm0 ); // O|Psi>
				NHAM.Multiply(1.0, Jnm0, 0.0 , Jnm1 );
				axpy(dim ,-bandCenter/HalfWidth,Jnm0, Jnm1 ); 
				for( int n = 0 ; n < numMoms1; n ++)
				{
					mu(m,n) += dot(dim ,Psi, Jnm0 );
					JT=Jnm0; Jnm0=Jnm1; Jnm1=JT;
					NHAM.Multiply(2.0, Jnm0, -1.0 , Jnm1 );   // J1 = 2*H*J0
					axpy(dim ,-2*bandCenter/HalfWidth,Jnm0, Jnm1 ); 
				}
				JT=Jm0; Jm0=Jm1; Jm1=JT;
				NHAM.Multiply(2.0, Jm0, -1.0 , Jm1 );   // J1 = 2*H*J0
				axpy(dim ,-2*bandCenter/HalfWidth,Jm0, Jm1 ); 
			}
		}
	}


	std::ofstream outputfile( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom0+"x"+sNumMom1+"RV"+snumRV+".chebmom2D").c_str() );
	typedef std::numeric_limits< double > dbl;
	outputfile.precision(dbl::digits10);
        outputfile<<dim<<" "<<2*HalfWidth<<"  "<<bandCenter<< std::endl;
        outputfile<<numMoms0<<" "<<numMoms0<< std::endl;
	for( int m = 0 ; m < numMoms0 ; m++)
	for( int n = 0 ; n < numMoms1 ; n++)
	{
		double scal = 4.0/numSV/numBlock;
		if(m == 0) scal *= 0.5;
		if(n == 0) scal *= 0.5;

		complex mu0 = scal*mu(m,n); 
		
		outputfile<<mu0.real()<<" "<<mu0.imag()<<std::endl;
	}
	outputfile.close();

	delete [] data;

	return 0;
}

