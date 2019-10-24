#ifndef QUANTUM_STATES
#define QUANTUM_STATES

#include <fstream>
#include <vector>
#include <complex>
#include <numeric>
using namespace std;  // for std::complex, and std::vector

struct StateParams
{

    int numStates;
    int dimension;
    const string sequence_type;
    const string origin_type;
}
class State
{
    public:
        void allocate(const int dim){
            data = vector<complex>(dim);
        };

        void release(){

        };
        int dimension(){ return data.size() };

    //Virtual methods
        virtual void create() = 0;
        virtual complex* initialPtr() = 0 ;

    private:
        vector<complex> data;
}

class RandomPhase: public State{

    virtual void create(){

    vector<complex>::iterator it; 
    for (it=data.begin();i1!=data.end();++it)
    {
        const double phi = 2*M_PI*rand()/RAND_MAX;
        *it = complex( cos(phi), sin(phi) )/sqrt(data.size());
    } 
        

        
    };
    virtual complex* initialPtr() = 0 ;
}

class Localized : public State{

    virtual void create() = 0;
    virtual complex* initialPtr() = 0 ;
}

class Momentum : public State{

    virtual void create() = 0;
    virtual complex* initialPtr() = 0 ;
}

class StateCreator{

    public:
    void setStateType(State *s, const StateParams *par )
    {
        _state.config( par );
    }

    bool setNewState()
    {
        if( _state.dimension() == 0 )
            _state.allocate( _stateParams-> dimension );
            _state.setNumStates( _stateParams-> dimension );

            
        bool final_state =  _state.create();
    
        if ( !final_state )
            _state.release();
    
       return !final_state; 
    }

    complex* getStatePtr()
    {
        return _state.initialPtr():
    }

    private:
        State *_state;
        StateParams *_stateParams;
}

/**
 * @brief	The Random Phase is class which create a state initialized with exp(i phi )
            where phi a random number. Is mostly used for stochastic trace approximation.
 */
class RandomPhase : public State
{

}

class LocalPosition :  public State
{

}

class LocalMomentum :  public State
{

}

#endif
