#ifndef HOPPING_LIST_UTIL 
#define HOPPING_LIST_UTIL

#include <array>
#include <string>
#include <string>
#include <sstream>
#include <ctype.h>


inline int 
index_aliasing(const array<int, 3>& index,const array<int, 3>& bound )
{
    return ( (index[2]+bound[2])%bound[2] * bound[1] + (index[1]+bound[1])%bound[1] ) * bound[0] + (index[0]+bound[0])%bound[0] ;
}


inline 
array<int,5> tag_to_indices(const string& tag)
{
	const int max_indices = 5;
    array<int,max_indices> indices;

    int count = 0 ;
    stringstream ss(tag);
    for( auto& ti : indices)
    if ( count < max_indices )
    {
        ss>>ti;
        count++;
	}
	
return indices; 
}


#endif
