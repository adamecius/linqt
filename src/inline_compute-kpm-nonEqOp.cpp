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
	if ( !(argc == 7 || argc == 8) )
	{
		chebyshev::printHelpMessage();
		return 0;
	}
	else
		chebyshev::printWelcomeMessage();
	
	const std::string
		LABEL = argv[1],
		S_OPR = argv[2],
		S_OPL = argv[3],
		S_NUM_MOM = argv[4];


	const int numMoms= atoi(argv[4]);
	chebyshev::Moments2D chebMoms(numMoms,numMoms); //load number of moments

	SparseMatrixType OP[3];
	OP[0].SetID("HAM");
	OP[1].SetID(S_OPR);
	OP[2].SetID(S_OPL);

	// Build the operators from Files
	SparseMatrixBuilder builder;
	for (int i = 0; i < 3; i++)
	{
		std::string input = "operators/" + LABEL + "." + OP[i].ID() + ".CSR";
		builder.setSparseMatrix(&OP[i]);
		builder.BuildOPFromCSRFile(input);
	};
	//CONFIGURE THE CHEBYSHEV MOMENTS
	chebMoms.SystemLabel(LABEL);
	chebMoms.BandWidth( atof(argv[5]) );
	chebMoms.BandCenter( atof(argv[6]) );
	chebMoms.SystemSize(OP[0].rank() );
	chebMoms.Print();


	//Define thes states youll use
	//Factory state_factory ;

	//Compute the chebyshev expansion table
	srand(time(0));
	int num_states = 1 ;
	if( argc == 8)	num_states = atoi(argv[7]);
	
	chebyshev::CorrelationExpansionMoments(num_states, OP[0], OP[1], OP[2], chebMoms, RANDOM_STATE );

	//Save the table in a file
	std::string outputfilename="NonEqOp"+S_OPR+"-"+S_OPL+LABEL+"KPM_M"+S_NUM_MOM+"x"+S_NUM_MOM+".chebmom2D";
	chebMoms.saveIn(outputfilename);

	std::cout<<"End of program"<<std::endl;
	return 0;
}

