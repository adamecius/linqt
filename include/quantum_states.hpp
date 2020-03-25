#ifndef QUANTUM_STATES
#define QUANTUM_STATES

#include <fstream>
#include <vector>
#include <complex>
#include <numeric>
using namespace std; // for std::complex<double> , and std::vector


enum StateType
{
	LOCAL_STATE=1,
	RANDOM_STATE=2
};

typedef std::complex<double>  Complex;
typedef std::vector<Complex>  Vector;

namespace qstates
{
	inline
	int FillWithRandomPhase(Vector& X)
	{
		for(auto& elem : X )
		{
			auto phi = 2.0*M_PI*(double)rand() / (double)RAND_MAX ;
			elem = Complex(cos(phi)/X.size(),sin(phi)/X.size() );
		}
		return 0;	
	}
}
#endif
