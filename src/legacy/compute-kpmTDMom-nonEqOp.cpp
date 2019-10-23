
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



int GetMomFromTol(const double tol, const double bandWidth, const double broadening )
{
	
	const int MAXMom = 100000;
	const double nb = 2.0*broadening/bandWidth/3.0;//The factor 3 is to take into account for 1/x to be at x= 3; 
	const double Xmin=-0.99,Xmax= -Xmin ;
	const int nx = (Xmax-Xmin)/(0.1*nb); 
	const double dX = (Xmax-Xmin ) / (nx-1);

	std::vector< complex > coeff(MAXMom,0);  
	std::vector< complex > aprox(nx,0);  
	
	for(int m = 0 ; m < MAXMom ; m++ )
	{
		coeff[m] = CPGF_Fun(m,0,nb );
		if( m == 0 ) coeff[0]*=0.5;
			
		double err = 0.0;
		for(  int n = 0 ; n < nx ; n++ )
		{
			const double x = Xmin + n*(Xmax-Xmin)/(nx-1);
			const complex exact = 1.0/complex(-x,nb);
			aprox[n] += coeff[m]*cos(m*acos(x) );
			err += sqrt( std::norm( ( aprox[n]-exact)/exact) ) ; 
		}
		err = err/nx;
		if ( err < tol )
			return m ;
	}
	return -1;
};



void NormEvolOpt(CSRMatrix& NHAM,  const int dim,  complex* aux, complex* X, const double delta_t )
{
	
	const complex I = complex(0.0,1.0);	
	complex *JT, *J0=&aux[0*dim], *J1=&aux[1*dim];


	copy(dim, X, J0 ); for(int i=0; i < dim ; i++) X[i]=0.0;
	NHAM.Multiply(1.0, J0, 0.0 , J1 );  // J1 = H*J0

	int n = 0;
	complex cheb_coeff  =  j0( delta_t );
	while ( sqrt( std::norm(cheb_coeff)) > 1E-10 )
	{
		axpy(dim,cheb_coeff,J0,X);  // Psi = FL0(E)J0
		n++;
		cheb_coeff  =2.*pow(-I, n)*jn( n, delta_t );
	//	std::cout<<delta_t<<" "<<cheb_coeff<<" "<<n<<" "<<std::endl; 
		JT=J0; J0=J1; J1=JT;
		NHAM.Multiply(2.0, J0, -1.0 , J1 );   // J1 = 2*H*J0
	}
}; 


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


	double bandWidth, bandCenter, broadening;
	cfg.lookupValue("BandWidth",bandWidth);
	cfg.lookupValue("BandCenter",bandCenter);
	cfg.lookupValue("Broadening",broadening); broadening=broadening/1000.;
	const double nb = broadening*2.0/bandWidth;

	double Emin, Emax, dE;
	cfg.lookupValue("Emin",Emin);
	cfg.lookupValue("Emax",Emax);
	cfg.lookupValue("dE",dE); 

	int numMoms, numRV, nt; double delta_t;
	cfg.lookupValue("TimeStep",delta_t);
	cfg.lookupValue("NumberOfTimeSteps",nt);
	cfg.lookupValue("NumberOfMoments",numMoms);
	cfg.lookupValue("NumberOfRandVec",numRV);
	std::stringstream ss;
	std::string sNumMom, snumRV,snt,sdelta_t;
	ss << numMoms; sNumMom = ss.str(); ss.str(std::string());
	ss << numRV; snumRV = ss.str(); ss.str(std::string());
	ss << delta_t; sdelta_t = ss.str(); ss.str(std::string());
	ss << nt; snt = ss.str(); ss.str(std::string());

	
//	numMoms =GetMomFromTol(0.1, bandWidth, broadening );// This computes the number of moments for getting a Green's function within 1% error. 
//	std::cout<<numMoms<<std::endl;
//	return 0;


	//READ NORMALIZED HAMILTONIAN
	CSRMatrix NHAM,OPL,OPR;
	NHAM.ReadCSRMatrix( "operators/"+LABEL+".HAM.CSR"); //Normalized Hamiltonian
	NHAM.Rescale( 2.0/bandWidth  );
	NHAM.Optimize(numMoms*numMoms); //OPTIMIZE MATRIX

	//READ OPERATOR LEFT
	OPL.ReadCSRMatrix( "operators/"+LABEL+"."+S_OPL+".CSR");
	OPL.Optimize(numMoms*numMoms);

	OPR.ReadCSRMatrix( "operators/"+LABEL+"."+S_OPR+".CSR");
	OPR.Optimize(numRV);
	


	//INITIALIZE KPM PARAMETERS
	const int numMoms0=numMoms,numMoms1=nt; 
	chebMom mu(numMoms0,numMoms1); 
	const int dim = NHAM.Dim();
	complex* data = new complex[ 6*dim + numMoms0*numMoms1 ];
	complex* JR  = &data[0*dim];
	complex* JL  = &data[1*dim];
	complex* JR0 = &data[2*dim];
	complex* JR1 = &data[3*dim];
	complex* aux = &data[4*dim];
	complex* JLV = &data[5*dim];
		  mu.pmu = &data[6*dim];
	complex* JT;

	//OUTPUT VARIABLE
	srand(1521);


	//TIME DEPENDENT PARAMETERS
	const double hbar = 1.0; //0.6582119624; // ev. fs
	const double wn   = 1.0;// bandWidth/2.0/hbar; 
	//INITIALIZE ITERATION
//	std::cout<<"Random "<<numRV<<std::endl;
	for( int r = 0 ; r < numRV ; r++) 
	{
		//DRAW A RANDOM VECTOR
		CreateRandomVector( dim, JL );
		OPR.Multiply(1.0, JL, 0.0 , JR ); 
		for( int n = 0 ; n < numMoms1; n ++)
		{
			//include velocity
			OPL.Multiply(1.0, JL, 0.0 , JLV ); // O|Psi>
			//CHEBYSHEV ITERATION
			copy(dim, JR, JR0 ); 
			NHAM.Multiply(1.0, JR0, 0.0 , JR1); 
			for( int m = 0 ; m < numMoms0; m ++)
			{
				mu(m,n) += dot( dim ,JLV, JR0 );
				JT=JR0; JR0=JR1; JR1=JT;
				NHAM.Multiply(2.0, JR0, -1.0 , JR1 );   // J1 = 2*H*J0
			}		
			NormEvolOpt(NHAM,dim,aux,JL,-wn*delta_t );  // U(t)|PsiL>
			NormEvolOpt(NHAM,dim,aux,JR,-wn*delta_t );  // U(t)|PsiL>
		}
	}

	std::ofstream outputfile( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom+"NTD"+snt+"RV"+snumRV+".TDmom2D").c_str() );
	for( int m = numMoms0-1 ; m >= 0 ; m--)
	for( int n = numMoms1-1 ; n >= 0 ; n--)
	{
//		complex mu0 = mu(m,n).real()*dim/numRV/M_PI/hbar*2.0/bandWidth; 
		complex mu0 = mu(m,n).real()*numRV/dim; 
		outputfile<<m<<" "<<n<<" "<<mu0.real()<<" "<<std::endl;
	}
	outputfile.close();



	outputfile.open( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom+"NTD"+snt+"RV"+snumRV+".ETA_COND").c_str() );
	for( double eta = nb*0.1; eta<=nb; eta+=0.1*nb )
	{
		double x = 0.01;
		double cond =0;
		for( int m = 0 ; m < numMoms0 ; m++)
		for( int n = 0 ; n < numMoms1 ; n++)
		{
			complex mu0 = mu(m,n).real()*dim/numRV/M_PI/hbar*2.0/bandWidth*2.0/bandWidth; 
			cond += -exp(-eta*wn*n*delta_t ) *delta_t* ChebWeigthR(m,x,eta).imag()* cos(m* acos(x) )*mu0.real() ; 
		}
		outputfile<<eta<<" "<<cond<<" "<<std::endl;
	}	
	outputfile.close();


	outputfile.open( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom+"NTD"+snt+"RV"+snumRV+".T_COND").c_str() );
	{
		double eta =nb;
		double x = 0.01;
		double cond =0;
		for( int n = 0 ; n < numMoms1 ; n++)
		{
			cond = 0.0;
			for( int m = 0 ; m < numMoms0 ; m++)
			{
				complex mu0 = mu(m,n).real()*dim/numRV/M_PI/hbar*2.0/bandWidth; 
				cond +=-exp(-eta*wn*n*delta_t ) *delta_t * ChebWeigthR(m,x,eta).imag()*mu0.real(); 
			}
			
			outputfile<<n*delta_t<<" "<<cond<<" "<<std::endl;
		}
	}
	outputfile.close();

	outputfile.open( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+sNumMom+"NTD"+snt+"RV"+snumRV+".COND").c_str() );
	for( double x = -0.9; x<=0.9; x=x+0.01 )
	{
		double eta =nb;
		double cond =0;
		for( int m = 0 ; m < numMoms0 ; m++)
		for( int n = 0 ; n < numMoms1 ; n++)
		{
			complex mu0 = mu(m,n).real()*dim/numRV/M_PI/hbar*2.0/bandWidth*2.0/bandWidth; 
			cond += -exp(-eta*wn*n*delta_t ) *delta_t * ChebWeigthR(m,x,eta).imag()*mu0.real(); 
		}

		outputfile<<x*bandWidth/2.0<<" "<<cond<<" "<<std::endl;
	}
	outputfile.close();


	delete [] data;
	return 0;
	
}
