
#include <iostream>
#include "irregular_hamiltonian.hpp"

///This function rescale the hamiltonian, between the , (-alpha,alpha)
void  IrregularHamiltonian::Rescale(	const my::real _Emin,const my::real _Emax, const my::real _cutoff )
{
	typedef my::SparseMatrix::InnerIterator SpMatIt;

	const my::real ScalE=2.*_cutoff/(_Emax-_Emin);
	const my::real ME=(_Emax+_Emin)/2;

	#pragma omp parallel for
	for (my::integer k=0; k<Hamiltonian_.outerSize(); ++k)
		for ( SpMatIt it(Hamiltonian_ ,k); it; ++it)
		{
			const my::integer
			i=it.col(),j=it.row();

			my::scalar
			v_ij= it.value();
			if(i==j)
				v_ij=v_ij-ME;

			v_ij=v_ij*ScalE;

			Hamiltonian_.coeffRef(i,j) = v_ij;                    // alternative: mat.coeffRef(i,j) += v_ij;
		}

};

///This function rescale the hamiltonian, between the , (-alpha,alpha)

void  IrregularHamiltonian::Multiply( 	const my::integer memSep, const my::integer dim, const my::real alpha,const my::scalar* x,
		const my::real beta , my::scalar* y)
{
//#pragma omp parallel  for
	  for(my::integer i=0;i<memSep*dim;i+=memSep)
		  y[i]=beta* y[i] ;

//#pragma omp parallel for
	for (my::integer k=0; k<Hamiltonian_.outerSize(); ++k)
		for ( my::SparseMatrix::InnerIterator it(Hamiltonian_ ,k); it; ++it)
		{
//			std::cout<<it.value()<<" "<<std::endl;
			y[it.row()*memSep]+= alpha*it.value()*x[it.col()*memSep];
		}
};


