
// Used for OPENMP functions
#include <omp.h>

// Parsing library
#include <libconfig.h++>



// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <sstream>
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/

// MKL LIBRARIES
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_spblas.h"
typedef std::complex<double> complex;
typedef MKL_INT integer;


//CUSTOM LIBRARIES
#include "kernel_functions.hpp"
#include "algebra_functions.hpp"


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


	double bandWidth, bandCenter;
	cfg.lookupValue("BandWidth",bandWidth);
	cfg.lookupValue("BandCenter",bandCenter);


	double E0, bmin, bmax, db;
	cfg.lookupValue("E0",E0);
	cfg.lookupValue("bmin",bmin); bmin=bmin/1000.;
	cfg.lookupValue("bmax",bmax); bmax=bmax/1000.;
	cfg.lookupValue("db",db);	  db = db/1000.;	

	int numMoms, numRV;
	cfg.lookupValue("NumberOfMoments",numMoms);
	cfg.lookupValue("NumberOfRandVec",numRV);
	std::stringstream ss;
	std::string sNumMom, snumRV, sE0;
	ss << numMoms; sNumMom = ss.str(); ss.str(std::string());
	ss << numRV; snumRV = ss.str(); ss.str(std::string());
	ss << E0; sE0 = ss.str();ss.str(std::string());


	

	//READ NORMALIZED HAMILTONIAN
	CSRMatrix NHAM,OPL,OPR;
	NHAM.ReadCSRMatrix( "operators/"+LABEL+".NHAM.CSR"); //Normalized Hamiltonian
	NHAM.Optimize(2*numMoms); //OPTIMIZE MATRIX

	//READ OPERATOR LEFT
	OPL.ReadCSRMatrix( "operators/"+LABEL+"."+S_OPL+".CSR");
	OPL.Optimize(2);

	OPR.ReadCSRMatrix( "operators/"+LABEL+"."+S_OPR+".CSR");
	OPR.Optimize(2);



	//INITIALIZE KPM PARAMETERS
	const int dim = NHAM.Dim();
	complex* data = new complex[ 5*dim ];
	complex* Psi= &data[0*dim];
	complex* JL = &data[1*dim];
	complex* JR = &data[2*dim];
	complex* J0 = &data[3*dim];
	complex* J1 = &data[4*dim];
	complex* JT;


	std::ofstream outputfile( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom+"E0"+sE0+"RV"+snumRV+"LambdaConv.dat").c_str() );

//	std::cout<<bmin<<" "<<bmax<<std::endl;
	for( double broadening = bmin; broadening < bmax ; broadening+=db )
	{
		//OUTPUT VARIABLE
		srand(1521);
		const
		double 
		nE0 = E0*2.0/bandWidth,
		nbr = broadening*2.0/bandWidth;
		double output_quantity = 0.0 ;
			
		//INITIALIZE ITERATION
		
		for( int r = 0 ; r < numRV ; r++) 
		{
			for(int i=0;i < 5*dim; i++)
				data[i] = 0.0;
				//DRAW A RANDOM VECTOR
				CreateRandomVector( dim, Psi );

			//CASE RIGHT SIDE
			{
				//CASE m=0;
				OPR.Multiply(1.0, Psi, 0.0 , J0 ); // J0 = OPL|Psi>
				axpy(dim,ChebWeigthR( 0, nE0,nbr),J0,JR);  // JL = FL0(E)J0

				//ALL OTHER m;
				NHAM.Multiply(1.0, J0, 0.0 , J1 );  // J1 = H*J0
				for( int m = 1 ; m < numMoms; m ++)
				{
					JT=J0; J0=J1; J1=JT;
					axpy(dim,ChebWeigthR(m , nE0, nbr),J0,JR);  // JL = FL0(E)J0
					NHAM.Multiply(2.0, J0, -1.0 , J1 );   // J1 = 2*H*J0
				}
			}

			//CASE LEFT SIDE SIDE
			{
				//CASE m=0;
				copy(dim, Psi, J0 );for(int i=0;i < dim; i++)	Psi[i] = 0.0;
				axpy(dim,ChebWeigthL(0,nE0, nbr),J0,Psi);  // Psi = FL0(E)J0

				//ALL OTHER m;
				NHAM.Multiply(1.0, J0, 0.0 , J1 );  // J1 = H*J0
				for( int m = 1 ; m < numMoms; m ++)
				{
					JT=J0; J0=J1; J1=JT;
					axpy(dim,ChebWeigthL( m, nE0,nbr),J0,Psi);  // Psi = FL0(E)J0
					NHAM.Multiply(2.0, J0, -1.0 , J1 );   // J1 = 2*H*J0
				}
				OPL.Multiply(1.0, Psi, 0.0 , JL ); // J0 = OPL|Psi>
			}
			output_quantity += dot(dim ,JL, JR ).imag();
		}
		outputfile<<broadening*1000<<" "<< output_quantity*4.0/bandWidth/bandWidth/numRV/(double)dim<<std::endl;
		std::cout<<broadening*1000<<" "<< output_quantity*4.0/bandWidth/bandWidth/numRV/(double)dim<<std::endl;
	}
	outputfile.close();
	delete [] data;
	return 0;
}




