#ifndef QUANTUM_STATES
#define QUANTUM_STATES

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <chrono>
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
	typedef std::complex<double> scalar;

	struct generator
	{
		public:

		generator():count(0),label("default"), num_states(1),kind(RANDOM_STATE){}

		int NumberOfStates() const
		{
			return num_states;
		}


		void NumberOfStates(const int n)
		{
			num_states = n;
			spos.reserve(n);
		}
		
		void SystemSize( const int n)
		{
			dim = n;
			out = std::vector< scalar >(dim);
		}

		bool getQuantumState()
		{
			if( count < num_states )
			{
				
				std::fill(out.begin(), out.end(), 0.0);

				switch (kind)
				{
					case LOCAL_STATE:
						out[ spos[count] ] = 1.0;
						break;
					case RANDOM_STATE:
						const double norm = sqrt(out.size());
						for( auto&x : out )
						{
							auto phi = 2.0*M_PI*(double)rand() / (double)RAND_MAX ;
							x = Complex(cos(phi)/norm,sin(phi)/norm );
						}
						break;
				}
				count++;
				return true;
			}			
			return false;
		}
		
		std::string StateLabel()
		{return label;
		}
	
		
		std::vector< scalar >& State()
		{
			return out;
		}
		

		std::string label;
		int count;
		std::vector< scalar > out;
		std::vector<int> spos;
		int dim;
		int num_states;
		StateType kind;
		
};
	
	
	
	
	
	inline 
	generator CreateStateSet(std::ifstream& file, std::string kind)
	{
		generator data;
		if( kind == "local")
		{
			data.kind = LOCAL_STATE;

			int ibuff;  
			file>>ibuff;
			data.NumberOfStates(ibuff);
			
			for ( int i = 0 ; i < data.num_states; i++ )
			{
				int ibuff;
				file>>ibuff;
				data.spos.push_back(ibuff);			
			}
			assert( data.spos.size() == data.num_states  );			
		}
			
	return data;
	};
	
	
	inline
	generator LoadStateFile( string filename)
	{
		std::ifstream infile(filename);
    
		std::string kind("random_phase");
		if( infile.good() )
			infile>>kind;
		else
			std::cout<<"State file not found, using default"<<std::endl;

		std::cout<<"Building "<<kind<<" states"<<std::endl;
		auto states = CreateStateSet(infile,kind);
		std::cout<<"Success"<<std::endl;
		states.label=filename;
		return states;
	}

	inline
	int FillWithRandomPhase(Vector& X)
	{
		int kpm_seed = time(0); 	
		if( getenv("KPM_SEED") ) 
			kpm_seed = std::stoi(string(getenv("KPM_SEED")));
		srand( kpm_seed );
		std::cout<<"Current seed is "<<kpm_seed<<std::endl;
		const double norm = sqrt(X.size());
		int i = 0;
		for(auto& elem : X )
		{
			auto phi = 2.0*M_PI*(double)rand() / (double)RAND_MAX ;
			elem = Complex(cos(phi)/norm,sin(phi)/norm );
			if( i < 10 )
				std::cout<<elem.real()<<" "<<elem.imag()<<std::endl;
			i++;
			
		}


		return 0;	
	}
}
#endif
