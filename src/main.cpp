// C & C++ libraries
#include <iostream> /* for std::cout mostly */
#include <string>   /* for std::string class */
#include <fstream>  /* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include "sparse_matrix.hpp"
#include "quantum_states.hpp"
#include "chebyshev_solver.hpp"

void printHelpMessage();

void printWelcomeMessage();

int main(int argc, char *argv[])
{
	if (argc != 7)
	{
		printHelpMessage();
		return 0;
	}
	else
		printWelcomeMessage();

	const std::string
		LABEL = argv[1],
		S_OPR = argv[2],
		S_OPL = argv[3],
		S_NUM_MOM = argv[4];

	chebyshev::Configure conf;
	conf.maxMemory = 70; //bites

	conf.TableSize.push_back(atoi(argv[4]));
	conf.TableSize.push_back(atoi(argv[4]));
	conf.scaleFactor = atof(argv[5]);
	conf.shift = atof(argv[6]);

	MKL_SparseType OP[3];
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

	//Define thes states youll use
	//Factory state_factory ;

	//Compute the chebyshev expansion table
	chebyshev::MomTable ctable(conf);
	chebyshev::CorrelationExpansionMoments(1, OP[0], OP[1], OP[2], ctable);
	/*

	//Save the table in a file
	std::ofstream outputfile( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+S_NUM_MOM+"x"+S_NUM_MOM+"RV1.chebmom2D").c_str() );
	ctable*/
	return 0;
}

void printHelpMessage()
{
	std::cout << "The program should be called with the following options: Label Op1 Op2 numMom scaleFactor shift" << std::endl
			  << std::endl;
	std::cout << "Label will be used to look for Label.Ham, Label.Op1 and Label.Op2" << std::endl;
	std::cout << "Op1 and Op2 will be used to located the sparse matrix file of two operators for the correlation" << std::endl;
	std::cout << "numMom will be used to set the number of moments in the chebyshev table" << std::endl;
	std::cout << "scaleFactor and shift will be used to rescale the hamiltonian from H to scaleFactor*H + shift" << std::endl;
};

void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute a table needed for expanding the correlation function in Chebyshev polynomialms" << std::endl;
};
