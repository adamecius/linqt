/*
 * hopping.hpp
 *
 *  Created on: Aug 20, 2016
 *      Author: jgarcia
 */

#ifndef HOPPING_HPP_
#define HOPPING_HPP_

//Defines the types in the program
#include "types_definitions.hpp"
#include <vector>

struct Hopping
{

	Hopping(const my::integer _final_idx,
			const my::integer _final_orb,
			const my::integer _final_spin,
			std::vector<my::integer> _Dr,
			my::scalar _val ):
			final_idx(_final_idx),
			final_orb(_final_orb),
			final_spin(_final_spin),
			Dr(_Dr), val(_val)
	{

	}

	my::integer  final_spin;
	my::integer  final_orb;
	my::integer  final_idx;
	std::vector<my::integer> Dr;
	my::scalar val;
};



#endif /* HOPPING_HPP_ */
