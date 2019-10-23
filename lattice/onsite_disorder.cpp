

#include "onsite_disorder.hpp"
#include <iostream>

namespace NumCal{
void AddAndersonDisorder(Lattice& _lat, const my::real W)
{
	const my::real epsilon=std::numeric_limits<my::real>::epsilon() ;
	if(W <= epsilon ) return ;

	typedef IrregularHamiltonian::Hopping Hopping;
	custom_random::uniform_real_dist RandomEnergy(-0.5*W, 0.5*W );

	int my_count=0;
	for(my::integer i0=0; i0 < _lat.CellsInDir(0)  ; i0++ )
		for(my::integer i1=0; i1 < _lat.CellsInDir(1)  ; i1++ )
			for(my::integer i2=0; i2 < _lat.CellsInDir(2)  ; i2++ )
				for(my::integer io=0; io < _lat.OrbitalNumber(); io++ )
				{
					my::real
					Ei=RandomEnergy(global_generator);
					for(my::integer is=0; is < _lat.SpinNumber()    ; is++ )
					{
						const my::integer
						k0=_lat.IndexesToIndex( i0, i1, i2, io, is);
						_lat.IrregHam().SetHopping( Hopping(k0,k0,Ei) );
					}
				}
};

void AddChargedPuddles(Lattice& _lat, const my::real p,const my::real U, const my::real Rc  )
{




	const my::real epsilon=std::numeric_limits<my::real>::epsilon() ;
	if(U <= epsilon ) return ;


	typedef IrregularHamiltonian::Hopping Hopping;		//The hopping which in this case is the triplet for the sparse matrix
	custom_random::uniform_real_dist prob_test(0.0,1.0);	//The random number used for probability test
	custom_random::uniform_real_dist RandomEnergy ( -0.5*U , 0.5*U );

	const my::integer
	estimated_impurities=( p*(2*_lat.TotalOfOrbitals()/_lat.TotalOfOrbitals() ) );//The two is to overstimate the value

	std::vector< std::vector<my::real> >
	imp_pos( estimated_impurities );						//The vector of impurity positions

	std::vector< my::real >
	this_imp_pos(2);										//a singl impurity position

	my::integer
	impurity_count= 0;

	// first we go through all the possible positions
	for(my::integer i2=0; i2 < _lat.CellsInDir(2)  ; i2++ )
		for(my::integer i1=0; i1 < _lat.CellsInDir(1)  ; i1++ )
			for(my::integer i0=0; i0 < _lat.CellsInDir(0)  ; i0++ )
				for(my::integer io=0; io < _lat.OrbitalNumber(); io++ )
				{
					if( prob_test(global_generator) < p )	//If the probability test is true
					{
						/** If the number of slots for impurities is smaller
						 *  than the current number of impurities in the system,
						 *  a push_back must be done. Else, only an assigment should be done
						 **/
						this_imp_pos[0]=A[0][0]*i0 + A[1][0]*i1 + io*Delta[0];
						this_imp_pos[1]=A[0][1]*i0 + A[1][1]*i1 + io*Delta[1];

						if( imp_pos.size() <= impurity_count )
							imp_pos.push_back(this_imp_pos);
						else
							imp_pos[impurity_count]=this_imp_pos;

						impurity_count=impurity_count+1;
					}
				}

	//  std::cout<<"The impurity count is "<<impurity_count<<" the expected value is: "<< estimated_impurities<<std::endl;

	if (imp_pos.size()!=impurity_count-1)	//Release the overstimated memory
		imp_pos.erase (imp_pos.begin()+ impurity_count,imp_pos.end());

	std::vector<my::scalar > temp_diagonal(_lat.TotalOfOrbitals(),0 );
	//We incorporate the onsite energy of the puddles.
	for(my::integer i2=0; i2 < _lat.CellsInDir(2)  ; i2++ )
		for(my::integer i1=0; i1 < _lat.CellsInDir(1)  ; i1++ )
			for(my::integer i0=0; i0 < _lat.CellsInDir(0)  ; i0++ )
				for(my::integer io=0; io < _lat.OrbitalNumber(); io++ )
					for(my::integer it0=-1; it0 <= 1; it0++ )	// tile index 1
						for(my::integer it1=-1; it1 <= 1; it1++ )	// tile index 2. This is used to take into account the periodic boundary condition
						{
							const my::integer I0=i0+_lat.CellsInDir(0)*it0 ;
							const my::integer I1=i1+_lat.CellsInDir(1)*it1 ;
							const my::real
							r[2]={
									A[0][0]*I0+A[1][0]*I1+io*Delta[0],
									A[0][1]*I0+A[1][1]*I1+io*Delta[1]
							};
							for(my::integer imp=0;imp< impurity_count; imp++)
							{
								const my::real
								dist= 	(r[0]-imp_pos[imp][0])*(r[0]-imp_pos[imp][0])+
								(r[1]-imp_pos[imp][1])*(r[1]-imp_pos[imp][1]);
								const my::real
								Ei=RandomEnergy(global_generator)*exp(-0.5*dist/Rc/Rc);

								//Created for periodic puddles
								if(Ei > epsilon )
									for(my::integer is=0; is < _lat.SpinNumber()    ; is++ )
									{
										//compute the site where the energy is going to be changed
										const int
										k0=_lat.IndexesToIndex( I0, I1, i2, io, is);
										temp_diagonal[k0]=temp_diagonal[k0]+Ei;
									}
							}
						}
	// std::cout<<"Setting the hamiltonian in the sparse matrix"<<std::endl;

	for(my::integer k0=0;k0<temp_diagonal.size();k0++)
		_lat.IrregHam().SetHopping( Hopping(k0,k0,temp_diagonal[k0]) );

};

}
