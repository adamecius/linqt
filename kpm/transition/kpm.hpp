/**\class KPM
 *
 * \brief A multi-purpose class for KPM calculations
 *
 * This class is aimed to served as a multi-purpose class for
 * KPM approximation of spectral quantities. It defines specialized
 * linear algebra operations and memory allocation used for an efficient
 * performance of large scale calculations.
 *
 * \note This class is still in development process, use it at your own risk!
 *
 * \version 1.0
 *
 */
#ifndef KPM_HPP
#define KPM_HPP


// Defines the types of the program
#include "types_definitions.hpp"
// Defines macros of parallelism
#include "kpm_parallel.hpp"

//Defines the headers for the lattice class
#include "lattice.hpp"

#include <cassert>
#include <vector>
#include <iostream>
#include <string>


namespace kpm
{
	static const my::real ONE=1,ZERO=0,TWO=1,NONE=-1;
	static const my::integer VECPERSET=3;
};

class ChebyshevSet
{
public:

	///The Main constructor
	/*!
	 * Take as parameters:\n
	 * The dimension of the vector space.\n\n
	 *
	 * It initialize the following parameters:\n\n
	 *
	 * Initialize the chebyshev vectors \n
	 * Initialize the truncation order in a default value 100 \n
	 * Initialize the minimal and maximal bounds in -1, 1 with a cutoff of 0.1 \n
	 * Initialize the energy bounds between -1 and 1 with a set up of 100 energy points
	 */
	ChebyshevSet(const std::string _setLabel, my::integer _vecDim ):
		setLabel_(_setLabel),
		vecDim_(_vecDim),
		trunc_order_(100),
		num_cheb_set_(2),
		kpm_linalg( KPMLinalg(_vecDim ,kpm::VECPERSET*NumOfSets() ) )
		{

		TotVec_= new my::scalar [NumOfSets()*3*Dimension() ];
		dE_= (Emax_-Emin_ )/NE_;
		scale_factor_= (max_bound_-min_bound_ )/( 2.*cutoff_ );

		//As default two set of chebyshev vectors are used
		ChebVecs_= std::vector< std::vector< my::scalar* > >(NumOfSets());
		for(my::integer i=0;i<NumOfSets();i++)
		{
			//Each of this sets has three chebyshev vectprs
			ChebVecs_[i]=std::vector< my::scalar* >(VECPERSET);
			//Which are pointing intercalated into a huge block of memory
			for(my::integer j=0;j<VECPERSET;j++)
				ChebVecs_[i][j]=&TotVec_[NumOfSets()*j+i];
		}
}


	~Kpm()
	{
		free(TotVec_);
	}


	/**************************SETTERS*************************************/
	void SetTruncOrder(const my::size_t _trunc_order);

	void SetEnergyScale(	const my::real _min_bound,
			const my::real _max_bound,
			const my::real _cutoff);

	void SetEnergyShift(	const my::real _min_bound,
			const my::real _max_bound,
			const my::real _cutoff);



	/**************************GETTERS*************************************/
	///Gets the dimension of the vector space
	inline
	my::integer Dimension() const
	{
		return dim_;
	}
	///Gets the total number of Chebyshev sets
	my::integer NumOfSets() const
	{
		return num_cheb_set_;
	}
	///Gets the truncation order of the system
	inline
	my::size_t TruncOrder() const
	{
		return trunc_order_;
	};
	/// Access a element of a chebyshev set
	inline
	my::scalar*&
	ChebVec(const my::integer i , const my::integer j )
	{
		return ChebVecs_[i][j];
	}
	/***********************Public Members**************************/
	///Converts the Normalized energy into energy
	my::real EnormToEn(const my::real Enorm) const;
	///Converts the energy into the normalized energy
	my::real EnToEnorm(const my::real En) const;
	///Computes the jackson kernel function
	my::real
	jackson_kernel(const my::real m);


	///The Chebyshev Iteration Function
	/*!
	 * This functions compute the iteration of the Chebyshev vectors
	 * following the rule:\n
	 * T_{n}(x)= 2x T_{n-1}(x) - T_{n-2}(x)
	 * It takes as parameters:\n
	 * The lattice vector used to extract the regular and irregular part of the Hamiltonian\n
	 * The set id of of the set  where the iteration will be performed
	 * The Chebyshev truncation order\n
	 * Initialize the truncation order in a default value 100 \n
	 */
	inline
	void RandomFill(const my::integer set_id, const my::integer ChebVecIDX)
	{
		const my::integer eff_dim=Dimension()*NumOfSets()*VECPERSET;
		const my::integer inc=VECPERSET*NumOfSets();

		for(my::integer i=0; i< eff_dim; i+=inc)
			ChebVec(set_id,ChebVecIDX)[i]=distribution(global_generator);
	}

	inline
	void PrintVector(const my::integer set_id, const my::integer ChebVecIDX)
	{
		const my::integer eff_dim=Dimension()*NumOfSets()*VECPERSET;
		const my::integer inc=VECPERSET*NumOfSets();

		for(my::integer i=0;i< eff_dim; i+=inc)
			std::cout<<i/NumOfSets()/VECPERSET<<" "<<ChebVec(set_id,ChebVecIDX)[i]<<std::endl;
	}

	///The Chebyshev Iteration Function
	/*!
	 * This functions compute the iteration of the Chebyshev vectors
	 * following the rule:\n
	 * T_{n}(x)= 2x T_{n-1}(x) - T_{n-2}(x)
	 * It takes as parameters:\n
	 * The lattice vector used to extract the regular and irregular part of the Hamiltonian\n
	 * The set id of of the set  where the iteration will be performed
	 * The Chebyshev truncation order\n
	 * Initialize the truncation order in a default value 100 \n
	 */
	inline
	void ChebyshevIteration(	/*Lattice& _lattice,*/const my::integer set_id, const my::integer n )
	{
		switch(n)
		{
		case 0:
		{

			const my::real one =1.0;
			const my::real zero=0.0;
	//		_lattice.RegHam().Multiply  (dim,shift, ONE, ChebVec(set_id,1),ZERO, ChebVec(set_id,2));
	//		_lattice.IrregHam().Multiply(dim,shift, ONE, ChebVec(set_id,1),ONE , ChebVec(set_id,2));
		} break ;
		default:
		{
			const my::real alpha= 2.0;
			const my::real beta =-1.0;
			const my::real one =1.0;

	//		_lattice.RegHam().Multiply  (dim,shift, TWO , ChebVec(set_id,2), NONE, ChebVec(set_id,1));
	//		_lattice.IrregHam().Multiply(dim,shift, TWO , ChebVec(set_id,2), ONE , ChebVec(set_id,1) );
			//Swap(set_id,1,2);
		} break;
		}

	}


	/**************************PRIVATE VARIABLES*************************************/


	void PrintKPMLoad()
	{
		std::cout	<<"The instance: \""<<LABEL<<"\" of the KPM class is using:"<<std::endl
					<<NumOfSets()<<" sets of 3 Chebyshev vectors"<<std::endl
					<<"and is using an approximated of "<< NumOfSets()*sizeof(TotVec_[0])*3*Dimension()/1048576<<" Mb of memory."<<std::endl;

	}
private:
	const std::string setLabel_;
	my::integer vecDim_, num_cheb_set_,trunc_order_;
	std::vector< std::vector< my::scalar* > >ChebVecs_;
	my::scalar* Chebt_;
	my::scalar* TotVec_;


};



#endif
