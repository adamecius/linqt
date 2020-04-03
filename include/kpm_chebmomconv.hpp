#ifndef KPM_CHEBMOMCONV_HELPER
#define KPM_CHEBMOMCONV_HELPER

// C & C++ libraries
#include <iostream> /* for std::cout mostly */

namespace chebyshev
{
	
inline
void printHelpMessage()
{
	std::cout << "The program should be called with the following options: Label Op1 Op2 numMom E0" << std::endl
			  << std::endl;
	std::cout << "Label will be used to look for Label.Ham, Label.Op1 and Label.Op2" << std::endl;
	std::cout << "Op1 and Op2 will be used to located the sparse matrix file of two operators for the correlation" << std::endl;
	std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
	std::cout << "E0 represents the energy in eV on which the convergence will be tested" << std::endl;
};

inline
void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute a correlation function in terms Chebyshev polynomials" << std::endl;
};


#endif
