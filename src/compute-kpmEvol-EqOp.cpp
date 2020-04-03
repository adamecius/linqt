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

int main(int argc, char *argv[])
{
	if ( !(argc == 8 || argc == 7 ) )
	{
		chebyshev::time_evolution::printHelpMessage();
		return 0;
	}
	else
		chebyshev::time_evolution::printWelcomeMessage();
	
	const std::string
		LABEL  = argv[1],
		S_OPL  = argv[2],
		S_OPR  = argv[3],
		S_NMOM = argv[4],
		S_NTIME= argv[5],
		S_TMAX = argv[6];

	const int numMoms = atoi(S_NMOM.c_str() );
	const int numTimes= atoi(S_NTIME.c_str() );
	const double tmax = stod( S_TMAX );

	chebyshev::MomentsTD chebMoms(numMoms, numTimes); //load number of moments

	SparseMatrixType OP[3];
	OP[0].SetID("HAM");
	OP[1].SetID(S_OPL);
	OP[2].SetID(S_OPR);

	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	for (int i = 0; i < 3; i++)
	{
		std::string input = "operators/" + LABEL + "." + OP[i].ID() + ".CSR";
		builder.setSparseMatrix(&OP[i]);
		builder.BuildOPFromCSRFile(input);
	
		if( i == 0 ) //is hamiltonian
		//Obtain automatically the energy bounds
		 spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
	};
	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebMoms.SystemLabel(LABEL);
	chebMoms.BandWidth ( (spectral_bounds[1]-spectral_bounds[0])*1.0);
	chebMoms.BandCenter( (spectral_bounds[1]+spectral_bounds[0])*0.5);
	chebMoms.TimeDiff( tmax/(numTimes-1) );
	chebMoms.SystemSize(OP[0].rank() );
	chebMoms.SetHamiltonian(OP[0]);
	chebMoms.Print();

	chebyshev::TimeDependentCorrelations( 1, OP[1], OP[2], chebMoms , RANDOM_STATE);

	
	std::string outputfilename="EvolEqOp"+S_OPL+"-"+S_OPR+LABEL+"KPM_M"+S_NMOM+".chebmomTD";	
	std::cout<<"Saving convergence data in "<<outputfilename<<std::endl;
	chebMoms.saveIn(outputfilename);
	std::cout<<"End of program"<<std::endl;
	return 0;
}

