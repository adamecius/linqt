
include "hopping_list_util.hpp"

int hop_utils::catchOrbitId( const string  orbId )
{
	if ( orbId.size() ==1 )
	{
		if( orbId.back() == 'x' )
			return all_x;
		if( orbId.back() == 'y' )
			return all_y;
	}
	bool is_number = true;
	for( auto x : orbId )
		is_number *= isdigit(x) ;
	
	if( is_number) 
		return stoi(orbId);
	
	return ignore;
		
};
