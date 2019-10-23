/*
 * irregular_hamiltonian.hpp
 *
 *  Created on: 04/08/2016
 *      Author: Jose Hugo Garcia Aguilar
 */

#ifndef IRREGULAR_HAMILTONIAN_HPP
#define IRREGULAR_HAMILTONIAN_HPP

#include "mpi_util.hpp"
///Defines the types of the program
#include "types_definitions.hpp"
///Defines the SparseMatrix class
#include "sparse_matrix.hpp"

///The Irregular Hamiltonian Class
/*! This class is the responsible of
 * store the irregular part of a Hamiltonian
 * in a sparse matrix format. It also contain
 * the implementation of the Hamiltonian*Vector product
 * which is extensively used in iterative programs
 * and also contain the implementation used to rescale the
 * spectrum between (-alpha,alpha) interval
 */
class IrregularHamiltonian
{
public:
      typedef sparse::Entry Hopping;
      ///The default constructor
      IrregularHamiltonian( ){}

      ///The initialization constructor  constructor
      IrregularHamiltonian( const my::integer ncol_,
			    const my::integer nrow_):
			      Hamiltonian_( my::SparseMatrix (ncol_ , nrow_))
      {
      }
      /********************PUBLIC METHODS******************************/
      inline void
      Reserve(const my::integer _size)
      {
	Hamiltonian_.Reserve(_size);
      }
      void
      Rescale( const my::real _Emin,const my::real _Emax, const my::real _cutoff );

      ///Set a hopping of the Hamiltonian
      void SetHopping( const Hopping _hop )
      {
	Hamiltonian_.AddNewRawEntry(_hop);
      }

      void
      Multiply(const my::integer memSep, const my::integer dim, const my::real alpha,const my::scalar* x, const my::real beta , my::scalar* y);

      void Refine()
      {
	Hamiltonian_.SetFromRawEntries();
      }

	my::SparseMatrix& Matrix()  
	{
		return Hamiltonian_;
	}

private:
      my::SparseMatrix Hamiltonian_;
};


#endif
