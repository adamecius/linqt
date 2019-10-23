/*
 * hopping_list.hpp
 *
 *  Created on: 26/08/2016
 *      Author: jgarcia
 */

#ifndef LATTICE_HOPPING_LIST_HPP_
#define LATTICE_HOPPING_LIST_HPP_

#include "mpi_util.hpp"
#include "types_definitions.hpp"
#include "hopping.hpp"

#include <iostream>
class HoppingList
{
public:
	int AddHopping(my::integer _init_orbit, Hopping _hop )
	{
		//Check if the new initial index already exists
		for( my::integer i=0; i< hop_indexes.size(); i++ )
			//if exits add the hopping into the same entry of the hop_indexes array
			// and finalize the call of the function
			if( hop_indexes[i] == _init_orbit )
			{
				hop_array  [i].push_back(_hop);
				return 0;
			}
		// if it doesn't exists before, add the index to the
		// hop_indexes array and create a slot for it in
		// the hop_array structures, the append the new hopping
		// into that structure and finalize the call
		hop_indexes.push_back(_init_orbit );
		hop_array.push_back( std::vector< Hopping >() );
		hop_array[ hop_indexes.size()-1].push_back(_hop );
		return 0;
	}

	void PrintHoppingList()
	{
		NumCal::cout<<NumCal::endl<<NumCal::endl;
		for( my::integer i=0;i< hop_indexes.size() ; i++)
		{
			NumCal::cout<<"\tInitial Orb Index: "<<hop_indexes[i]<<NumCal::endl;
			NumCal::cout<<"\tConnections:"<<NumCal::endl;
			for( my::integer h=0;h< hop_array[i].size() ; h++)
				NumCal::cout<<"\t\t"
						 <<hop_array[i][h].final_idx<<" "
						 <<hop_array[i][h].Dr[0]<<" "
						 <<hop_array[i][h].Dr[1]<<" "
						 <<hop_array[i][h].Dr[2]<<" "
						 <<hop_array[i][h].val<<NumCal::endl;
		}

	}

	my::integer
	InnerCellIndex(my::integer i)
	{
		return hop_indexes[i];
	}

	std::vector< Hopping >&
	HopOfOrbs(my::integer i)
	{
		return hop_array[i];
	}

	my::integer
	TotNumRegHops()
	{
		return hop_indexes.size();
	}

	void
	RescaleHopping(my::scalar scal)
	{
		for( my::integer i=0;i< hop_indexes.size() ; i++)
			for( my::integer h=0;h< hop_array[i].size() ; h++)
				hop_array[i][h].val=scal*hop_array[i][h].val;
	}

	void
	ShiftOnSite(my::scalar shift)
	{
		for( my::integer i=0;i< hop_indexes.size() ; i++)
			for( my::integer h=0;h< hop_array[i].size() ; h++)
				if(hop_indexes[i] == hop_array[i][h].final_idx )
					hop_array[i][h].val=hop_array[i][h].val - shift;
	}

private:
	std::vector< std::vector< Hopping > > hop_array;
	std::vector< my::integer > hop_indexes;
};



#endif /* LATTICE_HOPPING_LIST_HPP_ */
