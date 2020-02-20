#ifndef W2SP_ARG 
#define W2SP_ARG

#include <deque>
#include <array>
#include <string>
#include <iostream>
#include <cassert>
#include "operator_utils.hpp"
using namespace std;
//TODO: VALIDATE OPERATOR LIST
/*! This is a class which has as main role to handle input arguments of the program.
 *  It should: 
 *  (a) Read and validate the arguments from a command line. Done 
 *  (b) Read and validate the arguments from a config file. Done */
class W2SP_arguments
{
	
	public:
	//Read arguments from argc  and argv[]
	int ReadArguments( int argc, char* argv[] )
	{
		arguments = deque< string >(argv,argv+argc);
		program_name = arguments.front(); arguments.pop_front();  
		
		current_error="ERROR: The program: "+program_name+" should be called with arguments (LABEL, Dim0, Dim1, Dim2 ). Type help for more information. ";
		if( arguments.empty() ){ cout<<current_error<<endl; exit(-1);};
		label = arguments[0]; arguments.pop_front();  
	
		if( label == "help" )
			std::cout<<"wannier2sparse have to be called as:"<<endl
					 <<"wannier2sparse LABEL Dim0 Dim1 Dim2 [Operators]"<<endl
					 <<"LABEL is a string that defines the system label, and it will be used to look for the LABEL_hr.dat, LABEL.uc, LABEL.xyz and LABEL.dis"<<endl
					 <<"Dim0,Dim1,Dim2 are three integers defining the dimension of the supercell. Dimx>0."<<endl
					 <<"[Operators] is an optional parameter which defines the list of operators to be computed."<<endl
					 <<"If absent, only the hamiltonian will be computed."
					 <<"If 'all', all available operators will be computed"
					 <<"If a list of parameters, ie,  Vx Vy Vz. Only those will be computed"<<std::endl;
					 
			
		cellDim = array<int, 3>({1,1,1});
		for(int i = 0 ; i < 3 ; i++ )
		{
			assert( !arguments.empty() );
			cellDim[i] = stoi( arguments.front()); arguments.pop_front();  
		}
		//If no arguments. Assume user wants only hamiltonian
		//If arguments=all. Assume user wants all operators
		if( arguments.front() == "all" )
		{
			for( auto OP : oputil::AVAIL_OPS)
				operators.push_back(OP);
			arguments = deque< string >();
		}
		//If user submit a list. Pass this list
		while( !arguments.empty() )
		{
			operators.push_back( arguments.front() ); arguments.pop_front();
		}
		//TODO VALIDATE OPERATOR LIST


		return 0;
	}

		public: 
		array<int, 3> cellDim;
		string  label,program_name, current_error;
		deque< string > arguments; 	//deque is a structure which allows for pop_front()
		deque< string > operators; 
};

#endif
