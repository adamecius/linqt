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

void printValidatedInput(chebyshev::Configure& conf);

int main(int argc, char *argv[])
{
	if ( !(argc == 8) )
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
	conf.maxMemory = 1500; //bites

	conf.TableSize.push_back(atoi(argv[4])); //load number of moments
	conf.TableSize.push_back(atoi(argv[4])); //load number of moments
	conf.scaleFactor = atof(argv[5]);
	conf.shift = atof(argv[6]);

	printValidatedInput(conf);
	
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
	int 
	num_states = atoi(argv[7]);
	
	chebyshev::MomTable ctable(conf);
	chebyshev::LocalCorrelationExpansionMoments(num_states, OP[0], OP[1], OP[2], ctable);

	//Save the table in a file
	std::string outputfilename="NonEqOp"+S_OPR+"-"+S_OPL+LABEL+"KPM_M"+S_NUM_MOM+"x"+S_NUM_MOM+"RV1.chebmom2D";
	ctable.saveIn(outputfilename);

	std::cout<<"End of program"<<std::endl;
	return 0;
}

void printHelpMessage()
{
	std::cout << "The program should be called with the following options: Label Op1 Op2 numMom scaleFactor shift (optional) num_states" << std::endl
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


void printValidatedInput(chebyshev::Configure& conf)
{
       for(int i=0; i < conf.TableSize.size();i++)
	{
	    assert( conf.TableSize[i] >0 );
	    if( i==0) std::cout<<"We will use: ";
		std::cout<<conf.TableSize[i]<<" ";
	}
	std::cout<<"chebyshev moments."<<std::endl;
        assert( conf.scaleFactor > 0.0 );
	std::cout<<"The Chebyshev rescaling will be done using:"<<std::endl;
	std::cout<<"Scale Factor = " <<conf.scaleFactor<<std::endl;
	std::cout<<"Energy Shift = " <<conf.shift<<std::endl;
};
