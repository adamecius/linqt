//external libraries
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
// C libraries
// C++ libraries
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <complex>
#include <limits>
//CUSTOM LIBRARIES
#include "time_handler.hpp"  //For timeConst::INF

namespace besselCoeff
{
	typedef std::complex<double> complex;
	//This function allows for casting the exponential for very negative values of Re(z) by setting the function to zero
	complex safe_exp(const complex z );
	double  safe_exp(const double z );
	//This function allows for casting the bessel function for very small values of x by setting the function to zero
	double safe_jn(const int n,const  double x );
	//This functions returns the exact integration of the bessel coefficient when the integration goes to infinity.
	complex BesselCoeffTinf( const double m, const complex z );
	//This functions returns the exact integration of the bessel coefficient when the integration goes to infinity.
	complex DiffBesselCoeffTinf( const double m, const complex z );

	complex dGCoeff( const double m, const double x );

	complex GCoeff( const double m, const double x );


	//This is an auxiliary function which can be used to compute 
	// the integrals using qawo algorithm or to construct a function
	// for the real and imaginary parts entering qags algoritm
	double Fc(const double nt, void *p);
	//cos(imz*nt)*safe_exp(rez*nt)*
	//This is an auxiliary function which can be used to compute 
	// the integrals using qawo algorithm or to construct a function
	// for the real and imaginary parts entering qags algoritm
	double Fts(const double nt, void *p);

	//This is an auxiliary function which can be used to compute 
	// the integrals using qawo algorithm or to construct a function
	// for the real and imaginary parts entering qags algoritm
	double Ftc(const double nt, void *p);
	//cos(imz*nt)*safe_exp(rez*nt)*
	//This is an auxiliary function which can be used to compute 
	// the integrals using qawo algorithm or to construct a function
	// for the real and imaginary parts entering qags algoritm
	double Fs(const double nt, void *p);

	void BesselCArr(const double ntphi, const double ntmin, const double ntmax,  double nen , int numC,  complex* barr ) 
	{					//time, energy, and broadening are normalized	

		double neta = 1.0/ntphi ; if ( ntphi==timeConst::INF ) neta= 0.0;
		const complex I(0.0,1.0);
		const complex z = complex( nen, neta );
		const complex Iz = I*z;

		//We first initialize the integration workspace
		gsl_set_error_handler_off();
		const int LIMIT_SIZE = 10000;
		gsl_integration_workspace *workspace;
		workspace = gsl_integration_workspace_alloc(LIMIT_SIZE);
		double epsabs=1E-12,  epsrel=0.0; 	
		double error; complex result;
		int status;
		gsl_function gslF; 			//Initialize integration workspace		for (int  m = 0; m < 2 ; m++)
		for(int n = 0 ; n < numC; n++ )
		if( ntmax==timeConst::INF )
			barr[n] = BesselCoeffTinf( n, z);
		else
		{
			double params[3] = { n, Iz.real(), Iz.imag() };
			gslF.params = &params;
			//The most stable way to perform the integration is to divide
			//the Bessel function in intervals (t0,t1) , and this can be 
			//easily done by using its zeroes. We always start with t0=0 
			//and for the first interval we choose the first zero
			barr[n] = 0;
			int nzero = 1;
			double t0 = ntmin; 
			double t1 = gsl_sf_bessel_zero_Jnu(n,nzero);
			//Determine the first zero above t0
			while ( t1 <= t0 ) 
			{
				nzero++;
				t1 = gsl_sf_bessel_zero_Jnu(n,nzero);
			}
			//Then we start the iterative process
			while( t0 != ntmax )
			{
				if( t1 > ntmax ) t1 = ntmax; //If the high part of the interval is larger than tmax use tmax
				//The simulation parameters in GSL format
				//Compute real Part
				gslF.function = &Fc;
				status = gsl_integration_qag(&gslF,t0,t1, epsabs,epsrel,LIMIT_SIZE,GSL_INTEG_GAUSS61, workspace, &result.real(), &error);
				if( error> epsabs )
					std::cout<<"The real part of the integration did not converged"<<std::endl;
				//Compute imaginary Part
				gslF.function = &Fs;
				status = gsl_integration_qag(&gslF,t0,t1, epsabs,epsrel,LIMIT_SIZE,GSL_INTEG_GAUSS61, workspace, &result.imag(), &error);
				if( error> epsabs )
					std::cout<<"The imaginary part integration did not converged"<<std::endl;
				barr[n]+= result;
				nzero++;
				t0=t1; t1= gsl_sf_bessel_zero_Jnu(n,nzero);
			}
			barr[n]*=-2.0*I*pow(I,-n);
		}
		gsl_integration_workspace_free(workspace);

		barr[0]*=0.5;

	return ;    
	};


	void BesselBastinArr(const double ntphi, const double ntmin, const double ntmax,  double nen , int numC,  complex* barr ) 
	{					//time, energy, and broadening are normalized	

		double neta = 1.0/ntphi ; if ( ntphi==timeConst::INF ) neta= 0.0;
		const complex I(0.0,1.0);
		const complex z = complex( nen, neta );
		const complex Iz = I*z;

		//We first initialize the integration workspace
		gsl_set_error_handler_off();
		const int LIMIT_SIZE = 10000;
		gsl_integration_workspace *workspace;
		workspace = gsl_integration_workspace_alloc(LIMIT_SIZE);
		double epsabs=1E-12,  epsrel=0.0; 	
		double error; complex result;
		int status;
		gsl_function gslF; 			//Initialize integration workspace		for (int  m = 0; m < 2 ; m++)
		for(int n = 0 ; n < numC; n++ )
		if( ntmax==timeConst::INF )
			barr[n] = DiffBesselCoeffTinf( n, z);
		else
		{
			double params[3] = { n, Iz.real(), Iz.imag() };
			gslF.params = &params;
			//The most stable way to perform the integration is to divide
			//the Bessel function in intervals (t0,t1) , and this can be 
			//easily done by using its zeroes. We always start with t0=0 
			//and for the first interval we choose the first zero
			barr[n] = 0;
			int nzero = 1;
			double t0 = ntmin; 
			double t1 = gsl_sf_bessel_zero_Jnu(n,nzero);
			//Determine the first zero above t0
			while ( t1 <= t0 ) 
			{
				nzero++;
				t1 = gsl_sf_bessel_zero_Jnu(n,nzero);
			}
			//Then we start the iterative process
			while( t0 != ntmax )
			{
				if( t1 > ntmax ) t1 = ntmax; //If the high part of the interval is larger than tmax use tmax
				//The simulation parameters in GSL format
				//Compute real Part
				gslF.function = &Ftc;
				status = gsl_integration_qag(&gslF,t0,t1, epsabs,epsrel,LIMIT_SIZE,GSL_INTEG_GAUSS61, workspace, &result.real(), &error);
				if( error> epsabs )
					std::cout<<"The real part of the integration did not converged"<<std::endl;
				//Compute imaginary Part
				gslF.function = &Fts;
				status = gsl_integration_qag(&gslF,t0,t1, epsabs,epsrel,LIMIT_SIZE,GSL_INTEG_GAUSS61, workspace, &result.imag(), &error);
				if( error> epsabs )
					std::cout<<"The imaginary part integration did not converged"<<std::endl;
				barr[n]+= result;
				nzero++;
				t0=t1; t1= gsl_sf_bessel_zero_Jnu(n,nzero);
			}
			barr[n]*=-2.0*I*pow(I,-n);
		}
		gsl_integration_workspace_free(workspace);

		barr[0]*=0.5;

	return ;    
	};



double Fc(const double nt, void *p)
{
	const double 
	nj  =( (double *)p)[0],
	rez =( (double *)p)[1],
	imz =( (double *)p)[2];	
	//return the function
	return cos(imz*nt)*safe_exp(rez*nt)*safe_jn((int)nj,nt);  
};

double Fs(const double nt, void *p)
{
	const double 
	nj  =( (double *)p)[0],
	rez =( (double *)p)[1],
	imz =( (double *)p)[2];	
	//return the function
	return sin(imz*nt)*safe_exp(rez*nt)*safe_jn((int)nj,nt);  
};

double Ftc(const double nt, void *p)
{
	const double 
	nj  =( (double *)p)[0],
	rez =( (double *)p)[1],
	imz =( (double *)p)[2];	
	//return the function
	return cos(imz*nt)*safe_exp(rez*nt)*nt*safe_jn((int)nj,nt);  
};

double Fts(const double nt, void *p)
{
	const double 
	nj  =( (double *)p)[0],
	rez =( (double *)p)[1],
	imz =( (double *)p)[2];	
	//return the function
	return sin(imz*nt)*safe_exp(rez*nt)*nt*safe_jn((int)nj,nt);  
};


double safe_exp(const double x )
{
	gsl_sf_result exp_res;

	if( gsl_sf_exp_e(x,&exp_res) == GSL_EUNDRFLW ) 	//If underflow return 0
		return 0.0;	
	else
		return exp_res.val;  
};

complex safe_exp(const complex z )
{
	const double x = z.real();
	const double y = z.imag();
	gsl_sf_result exp_res;

	if( gsl_sf_exp_e(x,&exp_res) == GSL_EUNDRFLW ) 	//If underflow return 0
		return 0.0;	
	else
		return exp_res.val*complex( cos(y), sin(y) );  
};

double safe_jn(const int n,const  double x )
{
	int status;
	gsl_sf_result Jn_res;
	if( n == 0 )
		status = gsl_sf_bessel_J0_e(x,&Jn_res);
	else
		status = gsl_sf_bessel_Jn_e(n,x,&Jn_res);
	
	if ( status == GSL_EUNDRFLW ) 	//If underflow return 0
		return 0.0;	
	else
		return Jn_res.val;  
};

complex BesselCoeffTinf( const double m, const complex z )
{ 
	const complex I   = complex( 0, 1.0	  );
	const complex fz  = sqrt( 1.0 - z*z );
	return  -I*pow( z - I*fz, m )/fz ; 
};

complex DiffBesselCoeffTinf( const double m, const complex z )
{ 
	const complex I   = complex( 0, 1.0	  );
	const complex fz  = sqrt( 1.0 - z*z );
	return  -I*( z + I*m*fz  )*pow( z - I*fz, m )/fz/fz/fz ; 
};

complex GCoeff( const double m, const double x )
{
	const complex I   = complex( 0, 1.0       );
	const complex fx  = sqrt( 1.0 - x*x );
	return  -I*pow( x - I*fx, m )/fx ;

};

complex dGCoeff( const double m, const double x )
{ 
	const complex I   = complex( 0, 1.0	  );
	const complex fx  = sqrt( 1.0 - x*x );
	return  -I*( x + I*m*fx  )*pow( x - I*fx, m )/fx/fx/fx ; 

};


};


