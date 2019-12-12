#ifndef CHEB_MOM
#define CHEB_MOM

#include <array>
#include <vector>
#include <complex>
#include <string>
#include <cassert>
#include <iostream>		/* for std::cout mostly */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/

using namespace std;

namespace chebyshev 
{
class Moments2D
{
	public: 
	Moments2D():numMoms({0,0}){};

	Moments2D(const int m0,const int m1):numMoms({m0,m1}), mu( vector<complex<double> >(numMoms[1]*numMoms[0], 0.0) ){};

	Moments2D( std::string momfilename );

	//GETTERS
	inline
	int SystemSize() const { return system_size; };

	array<int, 2> MomentNumber() const { return numMoms; };

	int HighestMomentNumber(const int i) const { return numMoms[i]; };

	inline
	int HighestMomentNumber() const { return  (numMoms[1] > numMoms[0]) ? numMoms[1] : numMoms[0]; };

	inline
	string SystemLabel() const { return system_label; };

	inline
	double BandWidth() const { return band_width; };

	inline
	double HalfWidth() const { return BandWidth()/2.0; };

	inline
	double BandCenter() const { return band_center; };


	//SETTERS
	inline
	void SystemSize(const int dim)  { system_size = dim; };

	inline
	void SystemLabel(string label)  { system_label = label; };

	inline
	void BandWidth( const double x)  { band_width = x; };

	inline
	void BandCenter(const double x) { band_center = x; };

	//OPERATORS
	inline
	complex<double>& operator()(const int m0,const int m1) { return mu[ m0*numMoms[1] + m1 ]; };

	inline
	bool operator == (const Moments2D& other) const  
	{ 
		return true;
			this->system_label	== other.system_label &&
			this->system_size	== other.system_size &&
			this->numMoms  		== other.numMoms  &&
			this->band_width 	== other.band_width && 
			this->band_center	== other.band_center&&
			this->mu 			== other.mu;
	};


	//Transformation
	void ApplyJacksonKernel()
	{
		const double
		phi_J0 = M_PI/(double)(numMoms[0]+1.0),
		phi_J1 = M_PI/(double)(numMoms[1]+1.0);
		
		double g_D_m0,g_D_m1;
		for( int m0 = 0 ; m0 < numMoms[0] ; m0++)
		{
			g_D_m0=( (numMoms[0]-m0+1)*cos( phi_J0*m0 )+ sin(phi_J0*m0)*cos(phi_J0)/sin(phi_J0) )*phi_J0/M_PI;
			for( int m1 = 0 ; m1 < numMoms[1] ; m1++)
			{ 
				g_D_m1=( (numMoms[1]-m1+1)*cos( phi_J1*m1 )+ sin(phi_J1*m1)*cos(phi_J1)/sin(phi_J1) )*phi_J1/M_PI;
				this->operator()(m0,m1) *= g_D_m0*g_D_m1;
			}
		}
	};
	
	
	private:
	std::string system_label;
	int system_size;
	array<int, 2> numMoms;
	double band_width,band_center;
	std::vector< std::complex<double> > mu;

};

};
#endif 
