
#ifndef ONSITE_DISORDER_HPP
#define ONSITE_DISORDER_HPP

#include <vector>
#include "types_definitions.hpp"
#include "custom_random.hpp"
#include "lattice_geometry.h"
#include "lattice.hpp"

namespace NumCal{


void AddAndersonDisorder(Lattice& _lat, const my::real W);
void AddChargedPuddles(Lattice& _lat, const my::real p,const my::real U, const my::real Rc  );

};

#endif
