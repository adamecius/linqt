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
	if ( !(argc == 6 || argc == 7 ) )
	{
		chebyshev::convergence::printHelpMessage();
		return 0;
	}
	else
		chebyshev::convergence::printWelcomeMessage();
	
	const std::string
		LABEL  = argv[1],
		S_OP   = argv[2],
		S_NMOM = argv[3],
		S_NTIME= argv[4],
		S_TMAX = argv[5];

	const int numMoms = atoi(S_NMOM.c_str() );
	const int numTimes= atoi(S_NTIME.c_str() );
	const double tmax = stod( S_TMAX );

	chebyshev::MomentsTD chebMoms(numMoms, numTimes); //load number of moments

	SparseMatrixType OP[2];
	OP[0].SetID("HAM");
	OP[1].SetID(S_OP);

	// Build the operators from Files
	SparseMatrixBuilder builder;
	std::array<double,2> spectral_bounds;	
	for (int i = 0; i < 2; i++)
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
	chebMoms.SystemSize(OP[0].rank() );
	chebMoms.Print();

//	chebyshev::TimeEvolution( OP[0], OP[1], chebMoms);

	
	std::string outputfilename="EvolEqOp"+S_OP+LABEL+"KPM_M"+S_NMOM+".chebmomTD";	
	std::cout<<"Saving convergence data in "<<outputfilename<<std::endl;
/*
	std::ofstream outputfile( outputfilename.c_str() );
	for( int m =0; m < chebMoms.HighestMomentNumber(); m++ )
		outputfile << chebMoms(m).real() <<std::endl;
	outputfile.close();
	std::cout<<"End of program"<<std::endl;
*/
	return 0;
}

