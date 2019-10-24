// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <stdlib.h>
#include "sparse_matrix.hpp"
#include "chebyshev_solver.hpp"

template <typename T>
std::string to_string( T x ) 
{
	std::stringstream ss;
	ss << x;
	return ss.str();
}

void printHelpMessage();

void printWelcomeMessage();


int main(int argc, char *argv[])
{	
	if( argc != 8)
		{ printHelpMessage(); return 0; }
	else 
		printWelcomeMessage

	std::string
	LABEL = argv[1],
	S_OPR = argv[2],
	S_OPL = argv[3],
	S_NUM_MOM=arg[4];

	chebyshev::Configure conf;
	conf.TableSize.push_back( atoi(arg[4]) );
	conf.TableSize.push_back( atoi(arg[4]) );
	conf.scaleFactor  = atof(argv[5]);
	conf.shift        = atof(argv[6]);
	
	MKL_SparseType HAM;
	SparseMatrixBuilder builder;
	builder.setSparseMatrix(&HAM);
//	HAM.ReadCSRMatrix( "operators/"+LABEL+".HAM.CSR");
	builder.ConstructFromCOO(numRows,numCols,rows,cols,vals);

	MKL_SparseType OPR;
	SparseMatrixBuilder builder;
	builder.setSparseMatrix(&OPR);
//	HAM.ReadCSRMatrix( "operators/"+LABEL+".OPR.CSR");
	builder.ConstructFromCOO(numRows,numCols,rows,cols,vals);

	MKL_SparseType OPL;
	SparseMatrixBuilder builder;
	builder.setSparseMatrix(&OPL);
//	HAM.ReadCSRMatrix( "operators/"+LABEL+".OPL.CSR");
	builder.ConstructFromCOO(numRows,numCols,rows,cols,vals);


	//Define thes states youll use
	State::Creator stateCreator;

	//Compute the chebyshev expansion table
	chebyshev::MomTable ctable(conf);
	chebyshev::CorrelationExpansionMoments(&stateCreator,HAM,OPL,OPR,ctable);


	//Save the table in a file
	std::ofstream outputfile( ("NonEqOp"+S_OPR+LABEL+"KPM_M"+S_NUM_MOM+"x"+S_NUM_MOM+"RV1.chebmom2D").c_str() );
	ctable
	return 0;
}

void printHelpMessage()
{
	std::cout<<"The program should be called with the following options: Label Op1 Op2 numMom scaleFactor shift"<<std::endl;
	std::cout<<"Label will be used to look for Label.Ham, Label.Op1 and Label.Op2"<<std::endl;
	std::cout<<"Op1 and Op2 will be used to located the sparse matrix file of two operators for the correlation"<<std::endl;
	std::cout<<"numMom will be used to set the number of moments in the chebyshev table"<<std::endl;
	std::cout<<"scaleFactor and shift will be used to rescale the hamiltonian from H to scaleFactor*H + shift"<<std::endl;
};

void printWelcomeMessage()
{
	std::cout<<"This program will compute a table needed for expanding the correlation function in Chebyshev polynomialms"<<std::endl;
};
