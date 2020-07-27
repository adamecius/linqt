//!  Wannier2Sparse library. 
/*!
  This library map a Wannier Hamiltonian in the Wannier90 format
  into a supercell in real, momentum or eigenvector spaces which is
  stored as an sparse matrix.
*/
#ifndef WANNIER2SPARSE
#define WANNIER2SPARSE
 
#include <cstdio>
#include <string>
#include <array>
#include <deque>
#include "tbmodel.hpp"

//! A normal member taking two arguments and returning an integer value.
/*!
	\param The system label 
	\param An array of thee positive integer elements that define the supercell dimensions
	\param An vector consisting in operators. If empty, only the hamiltonian is generated. (See more information about other operators in )
	\return returns zero if the process was succesful, or -1 otherwise. 
	\sa QTstyle_Test(), ~QTstyle_Test(), testMeToo() and publicVar()
*/
int wannier2sparse(std::string label, array<int, 3> cellDim, std::deque< std::string > op_list);
 
#endif 
