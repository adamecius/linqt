// C & C++ libraries
#include <iostream> /* for std::cout mostly */
#include <string>   /* for std::string class */
#include <fstream>  /* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include <chrono>


#include "kpm_noneqop.hpp" //Message functions
#include "chebyshev_moments.hpp"
#include "sparse_matrix.hpp"
#include "quantum_states.hpp"
#include "chebyshev_solver.hpp"


namespace time_evolution
{
	void printHelpMessage();
	void printWelcomeMessage();
}


int main(int argc, char *argv[])
{
	if ( !(argc == 4 || argc == 5 ) )
	{
		time_evolution::printHelpMessage();
		return 0;
	}
	else
		time_evolution::printWelcomeMessage();
	
	const std::string
		LABEL  = argv[1],
		S_NTIME= argv[2],
		S_TMAX = argv[3];

	const int numTimes= atoi(S_NTIME.c_str() );
	const double tmax = stod(S_TMAX );

	chebyshev::MomentsTD chebMoms(1, numTimes); //load number of moments

	
	SparseMatrixType OP;
	OP.SetID("HAM");
	
	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	std::string input = "operators/" + LABEL + "." + OP.ID() + ".CSR";
	builder.setSparseMatrix(&OP);
	builder.BuildOPFromCSRFile(input);
	//Obtain automatically the energy bounds
	spectral_bounds = chebyshev::utility::SpectralBounds(OP);
	
	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebMoms.SystemLabel(LABEL);
	chebMoms.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0);
	chebMoms.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5);
	chebMoms.TimeDiff( tmax/(numTimes-1) );
	chebMoms.SetAndRescaleHamiltonian(OP);
	chebMoms.Print();

	//Compute the chebyshev expansion table
	
	chebyshev::TimeEvolutionTest(chebMoms);

	std::string outputName  ="timeEvolution_"+LABEL+".dat";

	std::cout<<"Saving the data in "<<outputName<<std::endl;

	std::ofstream outputfile( outputName.c_str() );
	for( int n = 0 ; n < chebMoms.MaxTimeStep(); n++)
	{
		outputfile << n * chebMoms.TimeDiff() <<" "<< chebMoms(0,n).real() << std::endl;
	}
	outputfile.close();

	std::cout<<"The program finished succesfully."<<std::endl;
	return 0;
}

void time_evolution::printHelpMessage()
	{
		std::cout << "The program should be called with the following options: Label numTimeSteps MaxTime num_states_file" << std::endl
				  << std::endl;
		std::cout << "Label will be used to look for Label.Ham" << std::endl;
		std::cout << "numTimeSteps  will be used to set the number of timesteps in the chebyshev table" << std::endl;
		std::cout << "TimeMax will be set the maximum time where the correlation will be evaluted " << std::endl;
	};

	inline
void time_evolution::printWelcomeMessage()
	{
		std::cout << "WELCOME: This program will compute a table needed for expanding the correlation function in Chebyshev polynomialms" << std::endl;
	};
