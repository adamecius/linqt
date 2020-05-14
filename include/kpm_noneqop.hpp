#ifndef KPM_NONEQ_HELPER
#define KPM_NONEQ_HELPER

// C & C++ libraries
#include <iostream> /* for std::cout mostly */

namespace chebyshev
{
	
inline
void printHelpMessage()
{
	std::cout << "The program should be called with the following options: Label Op1 Op2 numMom BandWidth BandCenter (optional) num_states" << std::endl
			  << std::endl;
	std::cout << "Label will be used to look for Label.Ham, Label.Op1 and Label.Op2" << std::endl;
	std::cout << "Op1 and Op2 will be used to located the sparse matrix file of two operators for the correlation" << std::endl;
	std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
	std::cout << "BandWidth and BandCenter will be used to rescale the hamiltonian from H to (2*H-BandCenter)/BandWidth" << std::endl;
};

inline
void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute a table needed for expanding the correlation function in Chebyshev polynomialms" << std::endl;
};

inline
bool GetBatchSize(int& batchSize)
{
	batchSize = 3;
	if(!getenv("BATCH_SIZE") || atoi(getenv("BATCH_SIZE")) <= 0) //evaluation left-to-right
		std::cout<<"\nEnviroment variable BATCH_SIZE not set or invalid.\nSet to your custom N value throught the command export BATCH_SIZE N. Using sequential"<<std::endl;
	else
	{
		batchSize = atoi(getenv("BATCH_SIZE"));
		return true;
	}	
	return false;
};

	namespace convergence
	{	
		inline
		void printHelpMessage()
		{
			std::cout << "The program should be called with the following options: Label Op1 Op2 numMom E0" << std::endl
					  << std::endl;
			std::cout << "Label will be used to look for Label.Ham, Label.Op1 and Label.Op2" << std::endl;
			std::cout << "Op1 and Op2 will be used to located the sparse matrix file of two operators for the correlation" << std::endl;
			std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
			std::cout << "Brodening(meV) represents the energy broadening in meV" << std::endl;
			std::cout << "E0 represents the energy in eV on which the convergence will be tested" << std::endl;
		};

		inline
		void printWelcomeMessage()
		{
			std::cout << "WELCOME: This program will compute a correlation function in terms Chebyshev polynomials" << std::endl;
		};

	};
};
#endif
