#ifndef QUANTUM_STATES
#define QUANTUM_STATES

#include <fstream>
#include <vector>
#include <complex>
#include <numeric>
using namespace std;  // for std::complex<double> , and std::vector

struct StateParams
{

    int numStates;
    int dimension;
    const string sequence_type;
    const string origin_type;
};

class State
{
};

class Factory{
    public:
    bool CreateNewState(){ return true; }
    vector<complex<double> > getNewState(){ return vector<complex<double> >();};
};

#endif
