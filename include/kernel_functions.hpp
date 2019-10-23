

#ifndef KERNEL_FUNCTION
#define KERNEL_FUNCTION
// KPM FUNCTIONS

double JacksonKernel(const double M, const double m0)
{
	const double theta = M_PI/( M + 1.0 );
	return ( ( M-m0+1.0)*cos(m0*theta) + sin(m0*theta)*cos(theta)/sin(theta) )/(M+1.0) ;
}


complex CPGF_Fun( const double m, const double x, double lambda )
{ 
	const complex I   = complex( 0, 1.0	  );
	const complex z   = complex( x, lambda);
	const complex fz  = sqrt( 1.0 - z*z );
	return  -2.0*I*pow( z - I*fz, m )/fz ; 
};

complex IncludeTemperature( complex* mom , const int numMom, const double mu, const double temp )
{
	
}


complex ChebWeigthR( const double m, const double x, double lambda= 1.0)
{
	const complex out = CPGF_Fun(m,x,lambda);
	if( m != 0 )
		return out;
	else
		return 0.5*out; 
};

complex ChebWeigthL( const double m, const double x, double lambda= 1.0)
{
	return ChebWeigthR(m,x,lambda);
/*	const double theta = m*acos(x);
	const complex out = CPGF_Fun(m,x,lambda); 
	if( m != 0 )
		return out;
	else
		return 0.5*out; 
*/		
};


complex BesselKernel( const int m, const double x, const double t, const double eta )
{
 	const double phi = x*t;
	const double alpha = exp(-eta*t)*jn( m, t );
	const complex I = complex(0,1); 
	const complex coeff =  complex( -sin(phi)*alpha, cos(phi)*alpha  );
	
	if( m == 0)
		return coeff*j0( t );
	else
	if ( m%2 == 0)
		if ( (m/2)%2 == 0)
			return  2.*coeff*jn( m, t );
		else
			return -2.*coeff*jn( m, t );

	else
		if ( ((m+1)/2)%2 == 0)
			return  2.*I*coeff*jn( m, t );
		else
			return -2.*I*coeff*jn( m, t );
};


double BesselKernel_Del( const int m, const double x, const double t, const double eta )
{
 	const double phi = x*t;
	const double alpha =-exp(-eta*t)*jn( m, t )/M_PI;
	
	if( m == 0)
		return cos(phi)*alpha;
	else
	if ( m%2 == 0)
		if ( (m/2)%2 == 0)
			return  2.*cos(phi)*alpha;
		else
			return -2.*cos(phi)*alpha;

	else
		if ( ((m+1)/2)%2 == 0)
			return -2.*sin(phi)*alpha;
		else
			return  2.*sin(phi)*alpha;
};

double BesselKernel_Del( const int m, const double x, const double eta, const double t , const double omega )
{
	const double hbar= 0.6582119624; // eV fs
 	const double tn  = t/hbar;
 	const double phi = x*tn;
	const double alpha =-exp(-eta*tn)*jn( m, omega*tn )/M_PI/hbar;
	
	if( m == 0)
		return cos(phi)*alpha;
	else
	if ( m%2 == 0)
		if ( (m/2)%2 == 0)
			return  2.*cos(phi)*alpha;
		else
			return -2.*cos(phi)*alpha;

	else
		if ( ((m+1)/2)%2 == 0)
			return -2.*sin(phi)*alpha;
		else
			return  2.*sin(phi)*alpha;
};

#endif
