#ifndef SPARSE_BLOCK_MATRIX_HPP
#define SPARSE_BLOCK_MATRIX_HPP
#include "types_definitions.hpp"
#include "dense_matrix.hpp"
#include <vector>
#include <complex>
#include <iostream>
#include <string.h>


namespace sparse
{
	template< typename T>
	class BlockMatrix 
	{
		public:
		BlockMatrix() {} ;

		BlockMatrix(const qt::dimension rank, const qt::dimension block_rank ) 
		{
			SetMatDim(rank,block_rank);
		}; 

	
		//PUBLIC INLINE FUNCTIONS
				
		inline qt::dimension
		BlockRank() const { return block_rank_; };

		inline qt::dimension 
		BlockNum() const { return block_num_; };
		
		inline qt::dimension 
		Rank() const { return rank_; };

		inline T* beginPtr() { return &data[0]; }
		inline T* endPtr() { return &data[data.size()-1]; }

		inline void
		SetMatDim(const qt::dimension rank, const qt::dimension block_rank )
		{ 
			block_rank_=block_rank;
			rank_=rank;
			block_num_  = rank/block_rank;
			data = std::vector<T>( block_num_*block_rank*block_rank );
		};

		//PUBLIC INLINE FUNCTORS 
		inline T 
		operator ()(const qt::index bId, const qt::index i, const qt::index j) const 
		{
			return data[ bId*BlockRank()*BlockRank() + i*BlockRank()+ j ];
		}
		inline T& 
		operator ()(const qt::index bId, const qt::index i, const qt::index j)  
		{
			return data[ bId*BlockRank()*BlockRank() + i*BlockRank()+ j ];
		}

		//PUBLIC FUNCTIONS
		void Multiply(	const qt::real a,const qt::complex* X,  
						const qt::real b, qt::complex* Y)
		const {
			for( qt::index bId=0; bId < BlockNum(); bId++)
			for( qt::index bi=0; bi < BlockRank(); bi++)
			{
				const qt::index 
				idx = ( bId*BlockRank() + bi )*BlockRank() , 
				i   =  bi + bId*BlockRank() ;					
				//  copy data to O for local aritmetic
				qt::complex Op[BlockRank()];
				memcpy ( Op, &data[idx], sizeof(qt::complex)*BlockRank() );

				qt::complex t= 0.0; // use for the intermediate multiplication stepts
				for( qt::index bj=0; bj < BlockRank(); bj++)
				{
					const qt::index 
					j   =  bj + bId*BlockRank() ;					
					
					// Sum_j Oij X_j
					const qt::complex x=X[j];
					t.real(t.real() + x.real()*Op[bj].real() - x.imag()*Op[bj].imag()) ;
					t.imag(t.imag() + x.real()*Op[bj].imag() + x.imag()*Op[bj].real()) ;
				}
				qt::complex y=Y[i];
				y.real( b*y.real() + a*t.real());
				y.imag( b*y.imag() + a*t.imag());
				Y[i]= y;
			 }
				
		}

		//PUBLIC FUNCTIONS
		void AddBlock(	const qt::index bId,const qt::dense_matrix<T>& Ab )
		{
			
			for( qt::index i=0; i < BlockRank(); i++)
			for( qt::index j=0; j < BlockRank(); j++)
			{
				const qt::index 
				idx = (bId*BlockRank() + i )*BlockRank() + j ;
				data[idx] = Ab(i,j);
			}	
		}

		void scale(qt::real a, qt::real b)
		{
			for( qt::index bId=0; bId < BlockNum(); bId++)
			for( qt::index i=0; i < BlockRank(); i++)
			for( qt::index j=0; j < BlockRank(); j++)
			{
				const myIndex idx = (bId*BlockRank() + i )*BlockRank() + j ; 
				if( i ==j )
					data[idx]-= b ;
				
				data[idx]*= a ;
			}
		}


		void Print() const
		{
			for( qt::index bId=0; bId < BlockNum(); bId++)
			{
				std::cout<<" The block: " << bId<<std::endl;
				std::cout<<std::endl;
				for( qt::index i=0; i < BlockRank(); i++)
				for( qt::index j=0; j < BlockRank(); j++)
				{
					const qt::index idx = (bId*BlockRank() + i )*BlockRank() + j ; 
					std::cout<<i+bId*BlockRank()<<" "<<j+bId*BlockRank()<<" "<<data[idx].real()<<" "<<data[idx].imag()<<std::endl;
					
				}
			}
		}
		
		public:
		std::vector<T> data;
		// This class assumes that the matrix 
		qt::dimension rank_,block_rank_,block_num_;
	};
}; 




#endif
