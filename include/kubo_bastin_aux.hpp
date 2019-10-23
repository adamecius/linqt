#ifndef KUBO_BASTIN_AUX
#define KUBO_BASTIN_AUX

#include <iostream>
#include <vector>
#include <complex>

#include <omp.h>

// GSL Libraries
#include <gsl/gsl_math.h> //
#include <gsl/gsl_integration.h>
#include <gsl/gsl_chebyshev.h> // Chebyshev Series
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_result.h>

//Discrete Cosine Transform
#include <fftw3.h>

namespace KuboBastin
{
	typedef std::complex<double> complex;


//AUXILIARY FUNCTIONS
//the fermi dirac distribution
double fermi_window( const double x )
{
	
//	if( fabs(x)> 40 )
//		return 0.0;
//	else
//	if( fabs(x)> 10 )
//		return -exp(-fabs(x));
//	else 
		return -exp(x)/( exp(x) + 1.0 )/( exp(x) + 1.0 );  	
};
	
//AUXILIARY FUNCTIONS
//the fermi dirac distribution
double fermi_dirac( const double x )
{
//	if( x< -40 )
//		return 1.0;
//	else
//	if( x> 10 )
//		return exp(-x);
//	else 
//	if( x > 200 )
//		return 0.0;
		return 1.0/( exp(x) + 1.0 );  	
};

//This functions samples the best point for performing
//a Chebyshev integration
double indexToEnergy( const int m , const int M )
{
	return cos( M_PI*( m + 0.5 )/(double)M ) ;
};

//This functions appears inside the kubo bastin formula
//when the exponential is expressed in chebyshev polynomials
//the phi is for allowing real and imaginary parts of  exp( I e t )
double weighted_j0( const double t, const double ener, const double eta, const double phi)
{
	return gsl_sf_bessel_J0(t)*exp(-t*eta)*cos(ener*t - phi);
}

//This functions appears inside the kubo bastin formula
//when the exponential is expressed in chebyshev polynomials
//the phi is for allowing real and imaginary parts of  exp( I e t )
double weighted_jn( const int n, const double t, const double ener, const double eta, const double phi)
{

	return gsl_sf_bessel_Jn(n,t)*exp(-t*eta)*cos(ener*t - phi);
}


//GSL FUNCTIONS

//This function appears in the kernel of the time integration
//in the Kubo-Bastin formula. This is the first where there is only 
//the zeroth-order bessel function. Higher order Bessel functions
//are going to be obtained recursively
//It has been generalized for
// an arbitrary power of t^power, because one may need to perform
//recursions where t^1, and t^0 appears.
double
GSL_F(double t, void *p);

//Compute the integration of the first coefficients
complex Gamma_n0_s(const double ener,const double eta, const double tmax, const double  power );


//The AlphaCoeff class allows for computing and calculating the
//The Alpha coefficients appearing in the Kubo Bastin formula. 

class AlphaCoeff
{
	public:

	AlphaCoeff(const int _mom0=1,const int _mom1=1):
		mom0(_mom0),mom1(_mom1),
		data(NULL)
	{
		//One need to allocate memory for two large arrays
		//needed to store the coefficients before
		//and after integration for speedup
		//and two addional arrays needed for recursive
		// methods
		data     = new complex[ 2*mom0*mom1 ];  
	};
	
	~AlphaCoeff()
	{
		if( data!= NULL )
			delete [] data;
	}

	complex& GammaCoeff( const int m0, const int m1 )
	{
		return data[ m1*mom0+m0];
	}

	complex& operator() ( const int m0, const int m1 )
	{
		return data[ (m1+mom1)*mom0+m0 ];
	}

	void ComputeGammaCoeff(const double tmax, const double eta );

	void ComputeGammaCoeff_FL(const double tmax, const double eta );

	void ComputeAlphaCoeff( const double mu, const double temp, const double dir = 1.0)
	{
		#pragma omp parallel 
		{
			fftw_r2r_kind kind = FFTW_REDFT10;
			double* dct_in  = fftw_alloc_real(2*mom0);
			double* dct_out = &dct_in[mom0];
			fftw_plan plan_dct = fftw_plan_r2r_1d(mom0,dct_in,dct_out,kind,FFTW_ESTIMATE);

			const int 
			tid = omp_get_thread_num(), 
			numThreads = omp_get_num_threads(),
			bsize= (mom1+numThreads-1 ) / numThreads;

			const int
			m1o = ( tid   ) * bsize,
			m1f = ( tid +1) * bsize;

			for(int m1=m1o ; m1< m1f; m1++ )
			if( m1 < mom1 )
			{
				for(int m0=0 ; m0< mom0; m0++ )
				{
					const double x  = dir*( indexToEnergy( m0 , mom0) - mu  )/temp;
					const double fd = fermi_dirac(x)/mom0; 
					dct_in[m0] = fd*GammaCoeff(m0,m1).real();
				}		
				fftw_execute(plan_dct); 

				for(int m0=0 ; m0< mom0; m0++ )
					if( m0 == 0 )
						(*this)( m0,m1).real( dir*0.5*dct_out[m0] );
					else
						(*this)( m0,m1).real( dir*1.0*dct_out[m0] );

				//Compute the imaginary  part
				for(int m0=0 ; m0< mom0; m0++ )
				{
					const double x  = dir*( indexToEnergy( m0 , mom0) - mu  )/temp;
					const double fd = fermi_dirac(x)/mom0; 
					dct_in[m0] = fd*GammaCoeff(m0,m1).imag();
				}
				fftw_execute(plan_dct); 
				for(int m0=0 ; m0< mom0; m0++ )
					if( m0 == 0 )
						(*this)(m0,m1).imag( dir*0.5*dct_out[m0]);
					else
						(*this)(m0,m1).imag( dir*1.0*dct_out[m0]);
			}
			fftw_free(dct_in);
			fftw_destroy_plan(plan_dct);
		}
	}

	void ComputeAlphaCoeff_FL( const double mu, const double temp)
	{
		#pragma omp parallel 
		{
			fftw_r2r_kind kind = FFTW_REDFT10;
			double* dct_in  = fftw_alloc_real(2*mom0);
			double* dct_out = &dct_in[mom0];
			fftw_plan plan_dct = fftw_plan_r2r_1d(mom0,dct_in,dct_out,kind,FFTW_ESTIMATE);

			const int 
			tid = omp_get_thread_num(), 
			numThreads = omp_get_num_threads(),
			bsize= (mom1+numThreads-1 ) / numThreads;

			const int
			m1o = ( tid   ) * bsize,
			m1f = ( tid +1) * bsize;

			for(int m1=m1o ; m1< m1f; m1++ )
			if( m1 < mom1 )
			{
				//Compute the real part
				for(int m0=0 ; m0< mom0; m0++ )
				{

					const double x  = ( indexToEnergy( m0 , mom0) - mu  )/temp;
					const double fd = fermi_window(x)/temp/mom0 ;
					dct_in[m0] = fd*GammaCoeff(m0,m1).real();
				}
				fftw_execute(plan_dct); 

				for(int m0=0 ; m0< mom0; m0++ )
					if( m0 == 0 )
						(*this)( m0,m1).real(0.5*dct_out[m0]);
					else
						(*this)( m0,m1).real(1.0*dct_out[m0]);

				//Compute the imaginary  part
				for(int m0=0 ; m0< mom0; m0++ )
				{
					const double x  = ( indexToEnergy( m0 , mom0) - mu  )/temp;
					const double fd = fermi_window(x)/temp/mom0 ;
					dct_in[m0] = fd*GammaCoeff(m0,m1).imag();
				}
				fftw_execute(plan_dct); 
				for(int m0=0 ; m0< mom0; m0++ )
					if( m0 == 0 )
						(*this)(m0,m1).imag(0.5*dct_out[m0]);
					else
						(*this)(m0,m1).imag(1.0*dct_out[m0]);
			}
			fftw_free(dct_in);
			fftw_destroy_plan(plan_dct);
		}
	};

	int numMom(){ return mom0; };
	int numTMom(){ return mom1; };
	
	private:
	complex *data;
	int mom0,mom1;


};



//----------------GSL FUNCTIONS AND METHODS ------------------------//

//This function appears in the kernel of the time integration
//in the Kubo-Bastin formula. This is the first where there is only 
//the zeroth-order bessel function. Higher order Bessel functions
//are going to be obtained recursively
//It has been generalized for
// an arbitrary power of t^power, because one may need to perform
//recursions where t^1, and t^0 appears.
double
GSL_F(double t, void *p) 
{
	const double 
	ener =  (double)( (double*)p )[0],//energy should be normalized
	eta  =  (double)( (double*)p )[1],//eta should be normalized
	power=  (double)( (double*)p )[2],//0,1,2,3,4,...
	phi  =  (double)( (double*)p )[3],// 0 for real part, M_PI_2 for imaginary part
	J0 	  = weighted_j0(t,ener,eta,phi);
	
	if( power == 0.0 )
		return J0;
	else 
	if (power == 1.0 )
		return t*J0;
	else 
		return pow(t,power)*J0;

};

//Compute the Gamma_n0^s coefficients. For more information 
// about what these coefficients are used for read the documentation
//they are essential for computing KuboBastin coefficients. 
complex Gamma_n0_s(const double ener,const double eta, const double tmax, const double  power )
{		

//	gsl_set_error_handler_off();						//The errors are going to be handled by mea
	double param[4] = { ener , eta , power ,  0.0 } ;	//The parameters need to be formared this way for GSL

	//Auxilary GSL variables
	int status = 0;
	const int max_steps = 1000;
	const double TOL = 1E-11;
	double result, error;//variables to store the integration output
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(max_steps);
    gsl_function F;
    F.function = &GSL_F;
    F.params = param;   
   
	//The complex integration is broken into real and imaginary parts
	// phi =0 real part, phi=M_PI_2 imaginary part.  
	complex output=complex(0,0);
    for(double phi = 0 ; ;) //First, the imaginary part is computed
    {
		//The integration of the Bessel function will be performed
		//between the zeroes of the Bessel function, starting from 
		//the one at t=0;
		double result= 0;
		double t0=0,t1=gsl_sf_bessel_zero_J0(1); 
		int n=1;
		while( t0 != tmax)
		{
			double partial_result= 0;
			param[3] =phi;
			status = gsl_integration_qag (	&F,   // function
											t0,   // from
											t1,   // to
											TOL,  // eps absolute
											0.0,  // eps relative
											max_steps,
											GSL_INTEG_GAUSS61,//Best accuracy for smooth functions
											w,
											&partial_result,
											&error);
											
			if (status != 0 )	//Check if there is any error in the integration
			{
				printf ("error: %s\n", gsl_strerror (status));
				std::cerr<<"While integrating between times:" <<t0<<","<<t1<<std::endl
						 <<"using params: "<<std::endl
						 <<"ener = "<<param[0]<<std::endl
						 <<"eta  = "<<param[1]<<std::endl
						 <<"power= "<<param[2]<<std::endl
						 <<"phi  = "<<param[3]<<std::endl;
				gsl_integration_workspace_free(w);
				return 0.0;
			}
			
			result+=partial_result;
			
		//	if ( std::fabs(result) < 1E-16 ) //If the integration result is very small break;
		//		t0 = tmax;	
			//If it is not small continue
			n++;
			t0 = t1;	
			t1 = gsl_sf_bessel_zero_J0(n);
			if( t1 > tmax )//If the new zero goes beyond the limits, set the limit as final value
				t1 = tmax;		
		}
		// Add the partial result to the real or imaginary part
		if( phi == 0 )
			output.real( result) ;
		if( phi == M_PI_2)
			output.imag( result) ;

		if ( phi == M_PI_2 ) break; //After the second iteration, this will break the for loop
		phi = M_PI_2; // Change from real to imaginary part. 
	}
	gsl_integration_workspace_free(w);
	return output;
};


void AlphaCoeff::ComputeGammaCoeff_FL(const double tmax, const double eta )
{
        if( tmax < 0.0 ) //The calculation is performed in finite time and the coefficients are known analytically
        {
                for ( int m0 = 0 ; m0 < mom0 ; m0++ )
                {
                        //It computes the coefficients for specifics
                        //energy distributed for an efficient chebyshev integration
                        const double ener = indexToEnergy( m0 , mom0);
                        //This variable includes the broadening and energy appearingin the 
                        //exponential
                        const complex z = complex(ener, eta );
                        //The coefficients are known analytically
                        const complex I(0,1);
                        const complex f1 = sqrt(1.0 - z*z );
                        const complex div_fact = 1.0/f1;
                        complex Gpow = 1.0;
                        for ( int m1 = 0 ; m1 < mom1 ; m1++ )
                        {
                                GammaCoeff(m0,m1) =-I*Gpow*div_fact;
                                Gpow *= z - I*f1;
                        }
                }
        return ;
        }
        std::cout<<"Not implemented"<<std::endl;
};
void AlphaCoeff::ComputeGammaCoeff(const double tmax, const double eta )
{
	if( tmax < 0.0 ) //The calculation is performed in finite time and the coefficients are known analytically
	{
		for ( int m0 = 0 ; m0 < mom0 ; m0++ )
		{
			//It computes the coefficients for specifics
			//energy distributed for an efficient chebyshev integration
			const double ener = indexToEnergy( m0 , mom0);
			//This variable includes the broadening and energy appearingin the 
			//exponential
			const complex z = complex(ener, eta );
			//The coefficients are known analytically
			const complex I(0,1);
			const complex f1 = sqrt(1.0 - z*z );
			const complex div_fact = 1.0/f1/f1/f1;
			complex Gpow = 1.0;
			for ( int m1 = 0 ; m1 < mom1 ; m1++ )
			{
				GammaCoeff(m0,m1) =-I*Gpow*( z + I*m1*f1 )*div_fact;
				Gpow *= z - I*f1;
			}
		}
	return ;
	}

	//The first loop is over the energy index m0
	#pragma omp parallel for
	for ( int m0 = 0 ; m0 < mom0 ; m0++ )
	{
		
		//These are auxiliary complex number that are going to be needed for 
		// to storage temporary energy dependent recursive vecotrs
		complex coeff0, coeff1;
		
		//It computes the coefficients for specifics
		//energy distributed for an efficient chebyshev integration
		const double ener = indexToEnergy( m0 , mom0);
		//This variable includes the broadening and energy appearingin the 
		//exponential
		const complex z = complex(-eta, ener);

		//The second loop is over the time index m1
		for ( int m1 = 0 ; m1 < mom1 ; m1++ )
		{
			//The first set of coefficients needs to be integrated numerically
			if( m1 == 0 )
				GammaCoeff(m0,m1) = Gamma_n0_s(ener,eta,tmax, 1.0 );

			else
			if( m1 == 1 )
			{
				//To get the second set of coefficients 
				//we first need to compute the coefficients of lower order
				//Gamma_n0_s(ener,eta,tmax, 0.0 )
				//Computed a weithed Bessel function with phi=0 for real and phi=M_PI_2 for imaginary
				complex wj0 = gsl_sf_bessel_J0(tmax)*exp(z*tmax);	
				//Then, auxiliary coefficients need to be computed 
				coeff0= Gamma_n0_s(ener,eta,tmax, 0.0 );
				coeff1= 1.0 - wj0  + z*coeff0;
				//then compute the real Gamma_n1_s(ener,eta,tmax, 0.0 )
				GammaCoeff(m0,m1) = -tmax*wj0 + coeff0 + z*GammaCoeff(m0,m1-1) ;
			}
			else //All other coefficients are computed recursively using the first two
			{
				//Computed a weithed Bessel function
				complex wjn =  gsl_sf_bessel_Jn(m1-1,tmax)*exp(z*tmax);
				//First compute the auxiliary bessel functons
				coeff0= coeff0- 2.0*( wjn - z*coeff1 );// coeff2 = coeff0 - 2*( weighted_jn() - z*coeff1 )
				complex swap = coeff0; coeff0 = coeff1; coeff1=swap;
				//then compute the real alpha coefficients
				GammaCoeff(m0,m1) = GammaCoeff(m0,m1-2)-2.0*( tmax*wjn - coeff0 - z*GammaCoeff(m0,m1-1) );
			}
		}		
	}

};



};


#endif
