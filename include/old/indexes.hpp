#ifndef TBLINALG_INDEXES_HPP
#define TBLINALG_INDEXES_HPP


//C libraries
#include <cassert>
//C++99 libraries
#include <vector>
//Custom Libraries
#include "types_definitions.hpp"



class Indexes
{
	public:

	//Class constructor 
	Indexes(const size_t _numIndexes):
	numIndexes_(_numIndexes), totalSize_(1)
	{
		idxArray_= std::vector<kpm::integer>( _numIndexes );
		dimArray_= std::vector<size_t>( _numIndexes );
		for( size_t i=0; i< _numIndexes ; i++ )
		{
			dimArray_[i]=1;
			idxArray_[i]=0;
		}
	}



	kpm::integer operator()( const  size_t pos) const
	{
		assert( pos < NumOfIndexes());
		return idxArray_[pos];
	}  

	kpm::integer& operator()( const  size_t pos)
	{
		assert( pos < NumOfIndexes());
		return idxArray_[pos];
	}  


	size_t Dim(const size_t pos) const 
	{
		assert( pos < NumOfIndexes());
		return dimArray_[pos];
	}

	size_t& Dim(const size_t pos) 
	{
		assert( pos < NumOfIndexes());			
		return dimArray_[pos];
	}

	size_t Size() const 
	{
		size_t totalSize=1;
		for( size_t i=0; i< NumOfIndexes() ; i++ )
		{
			totalSize*=Dim(i) ;
		}	
		return totalSize;
	}


	size_t NumOfIndexes() const
	{
		return numIndexes_ ;
	}
	
	kpm::integer ReduceToOne() const
	{
		kpm::integer singleIdx=idxArray_[0];
		
		for ( size_t pos=1 ; pos < NumOfIndexes(); pos++ ) 
			singleIdx = singleIdx*dimArray_[pos] + idxArray_[pos] ;

		return singleIdx;  
	}

	kpm::integer ScatterToAll(size_t idx ) 
	{
		for ( size_t pos=0 ; pos < NumOfIndexes() ; pos++ )
		{	const size_t inv_pos=NumOfIndexes()-1-pos;
			idxArray_[inv_pos] =( idx )%dimArray_[inv_pos];
			idx=idx/dimArray_[inv_pos];
		}
	}

	private:
		size_t numIndexes_;
		size_t totalSize_;
		std::vector<kpm::integer> idxArray_;
		std::vector<size_t> dimArray_;
};


#endif
