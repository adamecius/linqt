#ifndef KPM_BLOCK_MATRIX_HPP
#define KPM_BLOCK_MATRIX_HPP
#include "types_definitions.hpp"
#include <vector>
#include <complex>
#include <iostream>
#include <string.h>

typedef long int myint;
typedef unsigned long myIndex;

namespace kpm
{
	template< typename T>
	class BlockMatrix 
	{
		public:
		BlockMatrix() {} ;

		BlockMatrix(const myIndex rank, const myIndex block_rank ) 
		{
			SetMatDim(rank,block_rank);
		}; 

	
		//PUBLIC INLINE FUNCTIONS
				
		inline unsigned long 
		BlockRank() const { return block_rank_; };

		inline unsigned long 
		BlockNum() const { return block_num_; };
		
		inline unsigned long 
		Rank() const { return rank_; };

		inline T* beginPtr() { return &data[0]; }
		inline T* endPtr() { return &data[data.size()-1]; }

		inline void
		SetMatDim(const myIndex rank, const myIndex block_rank )
		{ 
			block_rank_=block_rank;
			rank_=rank;
			block_num_  = ( rank+block_rank-1 )/block_rank;
			data = std::vector<T>( block_num_*block_rank*block_rank );
		};

		//PUBLIC INLINE FUNCTORS 
		inline T 
		operator ()(const myIndex bId, const myIndex i, const myIndex j) const 
		{
			return data[ bId*BlockRank()*BlockRank() + i*BlockRank()+ j ];
		}
		inline T& 
		operator ()(const myIndex bId, const myIndex i, const myIndex j)  
		{
			return data[ bId*BlockRank()*BlockRank() + i*BlockRank()+ j ];
		}

		//PUBLIC FUNCTIONS
		void Multiply(	const kpm::real a,const kpm::complex* X,  
						const kpm::real b, kpm::complex* Y)
		const {
			for( myIndex bId=0; bId < BlockNum(); bId++)
			{
				kpm::complex t= 0.0; // use for the intermediate multiplication stepts
				for( myIndex i=0; i < BlockRank(); i++)
				{
					kpm::complex O[BlockRank()];

					const myIndex idx = (bId*BlockRank() + i )*BlockRank() ; 

					//  copy data to O for local aritmetic
					memcpy ( O, &data[idx], sizeof(O) );
					
					for( myIndex j=0; j < BlockRank(); j++)
					{
						// Sum_j Oij X_j
						const kpm::complex x=X[idx+j];
						t.real(t.real() + x.real()*O[j].real() - x.imag()*O[j].imag()) ;
						t.imag(t.imag() + x.real()*O[j].imag() + x.imag()*O[j].real()) ;
					}
					kpm::complex y=Y[idx];
					y.real( b*y.real() + a*t.real());
					y.imag( b*y.imag() + a*t.imag());
					Y[idx]= y;
				}
			}	
		}

		//PUBLIC FUNCTIONS
		void AddBlock(	const size_t bId,const kpm::dense_matrix<T>& Ab )
		{
			if ( Ab.Rank() == BlockRank() )
			for( size_t i=0; i < BlockRank(); i++)
			for( size_t j=0; j < BlockRank(); j++)
			{
				const myIndex 
				idx = (bId*BlockRank() + i )*BlockRank() + j ;
				data[idx] = Ab(i,j);
			}	
		}

		void scale(real a,real b)
		{
			for( myIndex bId=0; bId < BlockNum(); bId++)
			for( myIndex i=0; i < BlockRank(); i++)
			for( myIndex j=0; j < BlockRank(); j++)
			{
				const myIndex idx = (bId*BlockRank() + i )*BlockRank() + j ; 
				if( i ==j )
					data[idx]-= b ;
				
				data[idx]*= a ;
			}
		}

		void Print() const
		{
			for( myIndex bId=0; bId < BlockNum(); bId++)
			{
				std::cout<<" The block: " << bId<<std::endl;
				std::cout<<std::endl;
				for( myIndex i=0; i < BlockRank(); i++)
				for( myIndex j=0; j < BlockRank(); j++)
				{
					const myIndex idx = (bId*BlockRank() + i )*BlockRank() + j ; 
					std::cout<<i<<" "<<j<<" "<<data[idx].real()<<" "<<data[idx].imag()<<std::endl;
					
				}
			}
		}
		
		public:
		std::vector<T> data;
		// This class assumes that the matrix 
		unsigned long rank_,block_rank_,block_num_;
	};
}; 




#endif
