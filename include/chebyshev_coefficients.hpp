#ifndef  CHEBYSHEV_COEFFICIENTS
#define CHEBYSHEV_COEFFICIENTS

#include <complex>		/* for std::vector mostly class*/

using namespace std;

/*
 * 
void cdot(const int dim, const std::complex<double>* x,const std::complex<double>* y, std::complex<double>* result)
{
	*result = 0.0;
	for(int i=0; i < dim; i++ )
		*result = std::conj(x[i])*y[i];	
	return ;
}**/

inline
double delta_chebF(double x, const double m)
{
	const double f0 =  sqrt(1.0 - x*x ) ;
	const double fm =  cos(m*acos(x));
	return fm/f0/ M_PI;
};

inline 
complex<double> greenR_chebF(double x, const double m)
{
	const complex<double> I(0.0,1.0);
	const double f0 =  sqrt(1.0 - x*x ) ;
	const complex<double> fm =  pow( x - I*f0, m );
	return -I*fm/f0;
};

inline 
complex<double> DgreenR_chebF(double x, const double m)
{
	const complex<double> I(0.0,1.0);
	const double f0 =  sqrt(1.0 - x*x ) ;
	const complex<double> fm =  pow( x - I*f0 , m )*(x + I*m*f0);
	return -I*fm/f0/f0/f0;
};

/*
struct delta_Functor
{
	delta_Functor(double _theta): x(cos(_theta) ), f0( sqrt(1.0 - x*x) ) {};
    void operator()(double& m) const { m = cos(m*acos(x))*f0;  std::cout<<"I am f0" << f0<<std::endl; };
	const double x, f0;
};

struct greenR_Functor
{
	greenR_Functor(double _theta): x(cos(_theta) ), f0( sqrt(1.0 - cos(_theta)*cos(_theta) ) ),fm(1.0),m(0.0)  {};
	
    void operator()(std::complex<double>& y) 
    { 
		y = fm;
		this->fm *= (x - I * f0 );
		y   = fm/f0;
		m+=1.0;
	};
	
	void reset (){ m=0; };
	
	const double x, f0;
	double m;
	std::complex<double> fm;
};

*/
#endif 
