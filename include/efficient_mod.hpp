/*
 * efficient_mod.hpp
 *
 *  Created on: Aug 29, 2016
 *      Author: jgarcia
 */

#ifndef EFFICIENT_MOD_HPP_
#define EFFICIENT_MOD_HPP_

#include "types_definitions.hpp"

inline
my::integer EffMod(my::integer i, my::integer size)
{
	return (i+size )%size;
}

#endif /* EFFICIENT_MOD_HPP_ */
