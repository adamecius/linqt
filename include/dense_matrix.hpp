#ifndef DENSE_MATRIX_HPP
#define DENSE_MATRIX_HPP

#include "types_definitions.hpp"
#include <vector>
#include <complex>
#include <iostream>


namespace qt {
	
template< typename T>
class dense_matrix
{
	public:
	dense_matrix() {} ;

	dense_matrix(const qt::dimension _num_rows,const qt::dimension _num_cols )
	{
		num_rows_=_num_rows;
		num_cols_=_num_cols;
		data = std::vector<T>( _num_rows*_num_cols );
	};


	inline
	T operator ()(const qt::index i, const qt::index j) const
	{
		return data[ i*num_cols_ + j ];
	}

	inline
	T& operator ()(const qt::index i, const qt::index j)
	{
		return data[ i*NumOfCols() + j ];
	}

	void Print() const
	{

		for( qt::index i=0; i < NumOfRows() ; i++)
		{
			for( qt::index j=0; j < NumOfCols() ; j++)
			{
				std::cout<<data[ i*NumOfCols() + j ]<<"\t";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}

	inline
	void Transpose() 
	{
		std::vector<T> dataT = data;
		for( qt::index i=0; i < NumOfRows() ; i++)
		for( qt::index j=0; j < NumOfCols() ; j++)
			data[i*NumOfCols() + j ] =  dataT[ j*NumOfCols() + i ];
	}

	inline
	void HermitianConjugate() 
	{
		Transpose();
		for( qt::index i=0; i < data.size() ; i++)
			data[i] = std::conj(data[i]);
	}

	inline
	qt::dimension NumOfRows() const
	{
		return num_rows_;
	};

	inline
	qt::dimension NumOfCols() const
	{
		return num_cols_;
	};

	T* beginPtr()
	{
		return &data[0];
	}

	T* endPtr()
	{
		return &data[data.size()-1];
	}

	typedef typename std::vector<T>::iterator iterator;
	iterator begin() { return data.begin(); }
	iterator end()   { return data.end();   }

	public:
	std::vector<T> data;
	qt::dimension num_rows_,num_cols_;

};

}



#endif
