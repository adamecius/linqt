/*
 * This file is part of wannier2sparse.
 *
 * Developed as a tool of LinQT package.
 * 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#include <string>
#include <iostream>

#include "w2sp_arguments.hpp"
#include "wannier2sparse.hpp"
using namespace std;




/*! This is the main function for the executable  wannier2sparse */
int main( int argc, char* argv[]){

	W2SP_arguments  args;
	args.ReadArguments(argc,argv);

	cout<<"Using "<<args.label<<" as the system's identification label"<<endl
		<<"This label will be used to detect the label.xyz, label_hr.dat, label.stdis, and label.uc files"<<endl;

	wannier2sparse(args.label, args.cellDim, args.operators);
	
	/*
	tbmodel model;
	model.readOrbitalPositions(args.label+".xyz"); 	std::cout<<" finished"<<std::endl;
	model.readUnitCell(args.label+".uc"); 			std::cout<<" finished"<<std::endl;
	model.readWannierModel(args.label+"_hr.dat"); 	std::cout<<" finished"<<std::endl;
	model.readStaticDisorder(args.label+".stdis"); 	std::cout<<" finished"<<std::endl;


	std::cout<<"Creating Hamiltonian"<<std::endl;
	model.Hopping_List().wrap_in_supercell(args.cellDim).save_hopping_list_as_csr(args.label+".HAM.CSR");
	for( auto op : args.operators)
	{
		std::cout<<"Creating "<<op<<std::endl;
		model.WannierOperator(op).wrap_in_supercell(args.cellDim).save_hopping_list_as_csr(args.label+"."+op+".CSR");
	}	
	*/

	cout<<"The programa finished"<<std::endl;
return 0;}
