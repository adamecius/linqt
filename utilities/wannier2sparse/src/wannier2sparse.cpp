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


#include "wannier2sparse.hpp"
using namespace std;
int wannier2sparse(string label, array<int, 3> cellDim, deque< string > op_list)
{
	tbmodel model;
	model.readOrbitalPositions(label+".xyz");
	model.readUnitCell(label+".uc"); 		
	model.readWannierModel(label+"_hr.dat");
	model.readStaticDisorder(label+".stdis");

	//Create the Hamiltonian
	printf(" Creating the hamiltonian ... \n");
	model.Hopping_List().wrap_in_supercell(cellDim).save_hopping_list_as_csr(label+".HAM.CSR");
	for( auto op : op_list)
	{
		printf(" Creating the operator %s \n", op.c_str());
		model.WannierOperator(op).wrap_in_supercell(cellDim).save_hopping_list_as_csr(label+"."+op+".CSR");
	}	
	
	printf(" Wannier2Sparse finished succesfully \n");
return 0;}
