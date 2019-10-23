/**\class ChebsyhevSet
 *
 * \brief ChebyshevSet is an auxiliary class for the kpm principal class
 *
 * This class is to provide some important method for the containter defined as chebyshev set.
 * It defines how the set iterates itself, and a few operations for elements beloning to set
 * it also defines the set in an efficient way in memory.

 * \note This class is still in development process, use it at your own risk!
 *
 * \version 1.0
 *
 */
#ifndef CHEBYSHEV_SET_HPP
#define CHEBYSHEV_SET_HPP


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


namespace NumCal
{

/** \addtogroup kpm
 *  @{
 */

namespace kpm
{
static const my::real ONE=1,ZERO=0,TWO=2,NONE=-1;
static const my::integer VECPERSET=3;
};

class ChebyshevSet
{
public:

	///The Null Constructor
	ChebyshevSet():
	setLabel_("NullLabel"),
	vecDim_(1),
	num_cheb_set_(2)
	{
		memorySize_= 0;
		linearMemBlock_==NULL;
	}
	///The Main constructor
	/*!
	 * Take as parameters:\n
	 * The label of the set
	 * The dimension of the vector space.\n
	 */
	ChebyshevSet(const std::string _setLabel, my::integer _vecDim ):
		setLabel_(_setLabel),
		vecDim_(_vecDim),
		num_cheb_set_(2)
		{
			memorySize_= NumOfSets()*3*Dimension();
			///The linear memory blocks tries to allocate the appropiate ammount of memory
			linearMemBlock_=new my::scalar [ memorySize_ ];
			if(linearMemBlock_ == NULL)
			{
				std::cerr<<"Unable to allocate the ammount of requested memory: "
						 << sizeof(my::scalar)*memorySize_/1048576<<" Mb"<<std::endl
						 <<"the program will be terminated "<<std::endl;
				std::exit(-1);
			}
			//Initialize each of the requested chebyshev set, (by default 2)
			ChebVecs_= std::vector< std::vector< my::scalar* > >( NumOfSets() );
			for(my::integer i=0;i<NumOfSets();i++)
			{
			//Each of this sets has kpm::VECPERSET=3  Chebyshev vectors
			ChebVecs_[i]=std::vector< my::scalar* >(kpm::VECPERSET);
			//Which are setted intercalated into a huge block of memory
			for(my::integer j=0;j<kpm::VECPERSET;j++)
				ChebVecs_[i][j]=&linearMemBlock_[NumOfSets()*j+i];
			}
		}


	~ChebyshevSet()
	{
		free(linearMemBlock_);
	}


	/**************************SETTERS*************************************/
	void SetNumOfSets(const my::integer _set_num)
	{
		std::cout<<"This function is not yet implemented"<<std::endl;
	};


	/**************************GETTERS*************************************/
	///Gets the label of the set
	inline
	std::string Label() const
	{
		return setLabel_;
	}

	///Gets the dimension of the vector space
	inline
	my::integer Dimension() const
	{
		return vecDim_;
	}
	///Gets the total number of Chebyshev sets
	my::integer NumOfSets() const
	{
		return num_cheb_set_;
	}

	/// Returns the MemorySeparation for the sets
	inline
	my::integer
	MemSep() const
	{
		return NumOfSets()*kpm::VECPERSET;
	}

	/// Returns the jth vector of the ith  set
	inline
	my::scalar*&
	ChebVec(const my::integer i , const my::integer j )
	{
		return ChebVecs_[i][j];
	}
	/***********************Public Members**************************/
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
	void ChebyshevIteration( Lattice& lat, const my::integer set_id, const my::integer n )
	{
		switch(n)
		{
		case 0:
		{
			lat.ApplyHamiltonian(MemSep(), kpm::ONE, ChebVec(set_id,1),kpm::ZERO, ChebVec(set_id,2));
		} break ;
		default:
		{
			lat.ApplyHamiltonian(MemSep(), kpm::TWO, ChebVec(set_id,2),kpm::NONE, ChebVec(set_id,1));
			Swap(set_id,1,2);
		} break;
		}

	}

	inline
	void Swap(my::integer set_id,my::integer vec_id0, my::integer vec_id1)
	{
		Chebt_= ChebVec(set_id,vec_id0);
		ChebVec(set_id,vec_id0) = ChebVec(set_id,vec_id1);
		ChebVec(set_id,vec_id1) = Chebt_;
	}
	/**************************PRIVATE VARIABLES*************************************/
	inline
	void PrintVector(const my::integer set_id, const my::integer ChebVecIDX)
	{
		const my::integer eff_dim=Dimension()*NumOfSets()*kpm::VECPERSET;
		const my::integer inc=kpm::VECPERSET*NumOfSets();

		for(my::integer i=0;i< eff_dim; i+=inc)
			std::cout<<i/NumOfSets()/kpm::VECPERSET<<" "<<ChebVec(set_id,ChebVecIDX)[i]<<std::endl;
	}

	void PrintSetLoad()
	{
		std::cout	<<"The instance: \""<<setLabel_<<"\" of the Chebsyshev vector  class is using:"<<std::endl
				<<NumOfSets()<<" sets of 3 Chebyshev vectors"<<std::endl
				<<"and is using an approximated of "<< sizeof(linearMemBlock_[0])*memorySize_/1048576<<" Mb of memory."<<std::endl;

	}
private:
	const std::string setLabel_;
	my::integer vecDim_, num_cheb_set_,trunc_order_;
	std::vector< std::vector< my::scalar* > >ChebVecs_;
	my::scalar* Chebt_;
	my::scalar* linearMemBlock_;
	my::size_t  memorySize_;


};

/** @}*/

}

#endif
