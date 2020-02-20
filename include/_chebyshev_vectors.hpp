#ifndef CHEB_VEC
#define CHEB_VEC

#include <array>
#include <vector>
#include <complex>
#include <string>
#include <cassert>
#include <iostream>		/* for std::cout mostly */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/

using namespace std;

namespace chebyshev 
{
class VectorList
{
	public:
	typedef std::complex<double> value_t ;
	typedef std::vector<value_t> vector_t ;
	typedef std::vector<vector_t> vectorList_t ;

	VectorList(const int vecSize = 0 , const int numVec = 0):
	numVec_(vecSize),
	vecSize_(numVec),
	scale_(0.01),shift_(0.0),
	_J0(vector_t(vecSize)),_J1(vector_t(vecSize)),
	vecList_( vectorList_t(numVec) )
	{
		for(int m=0; m < numVec; m++)
			vecList_[m] = vector_t(vecSize);
	};

	//GETTERS
	inline
	int VectorSize()   const  { return vecSize_; }

	inline
	int NumOfVectors() const  { return numVec_; }
	

	inline
	vector_t&  VecCheb0()
	{
		return _J0; 
	};

	inline
	value_t&  VecCheb0(const int m)
	{
		return _J0[m]; 
	};

	inline
	vector_t&  VecCheb1()
	{
		return _J1; 
	};



	inline
	value_t&  VecCheb1(const int m)
	{
		return _J1[m]; 
	};

	inline
	vector_t&  GetVector(const  int m)
	{
		assert( m> 0 && m < VectorSize());
		return vecList_[ m ]; 
	};


	value_t&  GetVectorElem(const  int m,const  int n)
	{
		assert( n>0 && n <VectorSize() );
		return GetVector(m)[n]; 
	};

	inline
	double Shift() const { return shift_; }

	inline
	double Scale() const { return scale_; }

	//SETTERS
	void SetInitVectors( vector_t &J0, vector_t &J1 ){};
	

	inline 
	int SetShift( const double shift) { shift_ = shift; return 0;}

	inline 
	int SetScale( const double scale) { assert(scale> 0); scale_ = scale; return 0;}

	inline
	int  SetVector(const int m, const  vector_t  X )
	{
		assert( X.size() == VectorSize());
		GetVector(m)=X;
		return 0;
	};


	//heavier functions
	void SetInitVectors( SparseMatrixType &NHAM,const Moments::vector_t& T0 );

	void SetInitVectors( SparseMatrixType &NHAM, SparseMatrixType &OP ,const Moments::vector_t& T0 );

	int IterateAll( SparseMatrixType &NHAM );

	int Multiply( SparseMatrixType &OP );

	double MemoryConsumptionInGB();
	
	

	private:
	int numVec_;
	int vecSize_;
	double scale_,shift_;
	vectorList_t vecList_;
	vector_t _J0,_J1;
};
};
#endif 
