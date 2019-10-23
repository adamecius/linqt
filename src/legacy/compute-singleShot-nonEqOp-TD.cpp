
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



struct CSRMatrix
{

	void Multiply(const complex a,const complex* x,const complex  b, complex* y );

	integer Dim(){ return dim_;};

	bool ReadCSRMatrix( const std::string input);

	void Optimize( const int numMul )
	{
		mkl_sparse_set_mv_hint(A,SPARSE_OPERATION_NON_TRANSPOSE,descr, (MKL_INT)numMul );
		mkl_sparse_set_dotmv_hint(A,SPARSE_OPERATION_NON_TRANSPOSE,descr, (MKL_INT)numMul );
		mkl_sparse_optimize (A);
	}

	private:
		struct matrix_descr descr;
		sparse_matrix_t A;
//		char descr[6];		
		integer dim_;
		std::vector<complex> values;
		std::vector<integer> columns, rowIndex;
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


	double bandWidth, bandCenter, broadening;
	cfg.lookupValue("BandWidth",bandWidth);
	cfg.lookupValue("BandCenter",bandCenter);
	cfg.lookupValue("Broadening",broadening); broadening=broadening/1000.;


	double Emin, Emax, dE;
	cfg.lookupValue("Emin",Emin);
	cfg.lookupValue("Emax",Emax);
	cfg.lookupValue("dE",dE); 

	int numMoms, numRV;
	cfg.lookupValue("NumberOfMoments",numMoms);
	cfg.lookupValue("NumberOfRandVec",numRV);
	std::stringstream ss;
	std::string sNumMom, snumRV;
	ss << numMoms; sNumMom = ss.str(); ss.str(std::string());
	ss << numRV; snumRV = ss.str();

	

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


	//OUTPUT VARIABLE
	std::ofstream outputfile( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom+"RV"+snumRV+".dat").c_str() );
	for( double E0=Emin; E0 <= Emax; E0+=dE ){ srand(1521);
	
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
	//		std::cout<<r<<" "<< output_quantity/(double)(r+1)/(double)dim<<std::endl;
		}
		output_quantity *= 1.0/(double)numRV/(double)dim;
		outputfile<<E0<<" "<< output_quantity*4.0/bandWidth/bandWidth<<std::endl;
	}
	outputfile.close();
	delete [] data;
	return 0;
}

