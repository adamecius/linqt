#include "chebyshev_moments.hpp"


void chebyshev::Vectors::SetInitVectors( SparseMatrixType &NHAM,const Moments::vector_t& T0 )
{
	assert( NHAM.rank() == this->SystemSize()&& T0.size() == this->SystemSize() );
	ChebV0 = T0; ChebV1 = T0;
	NHAM.Multiply( ChebV0, ChebV1 );
};


void chebyshev::Vectors::SetInitVectors( SparseMatrixType &NHAM, SparseMatrixType &OP ,const Moments::vector_t& T0 )
{
	assert( OP.rank() == NHAM.rank() && NHAM.rank() == this->SystemSize()&& T0.size() == this->SystemSize() );

	ChebV0 = T0; ChebV1 = T0;
	OP.Multiply( ChebV1, ChebV0 );
	NHAM.Multiply( ChebV0, ChebV1 );
};


int chebyshev::Vectors::IterateAll( SparseMatrixType &NHAM )
{
	const int  dim = NHAM.rank();
	assert( dim == this->SystemSize() );
	
	linalg::copy(dim, &ChebV0[0],&this->operator()(0) );
	for(int m=1; m < this->HighestMomentNumber(); m++ )
	{
		linalg::copy(dim, &ChebV1[0],&this->operator()(m) );
		NHAM.Multiply(2.0,ChebV1,-1.0,ChebV0);
		ChebV0.swap(ChebV1);
	}
	return 0;
};

int chebyshev::Vectors::Multiply( SparseMatrixType &OP )
{
	assert( OP.rank() == this->SystemSize() );
	const int  dim = OP.rank();
		
	Moments::vector_t OPV( dim );
	for(int m=0; m < this->HighestMomentNumber(); m++ )
	{
		linalg::copy(dim, &this->operator()(m), &OPV[0] ); 
		OP.Multiply(&OPV[0], &this->operator()(m) );
	}
	return 0;
};

double chebyshev::Vectors::MemoryConsumptionInGB()
{
	return sizeof(value_t)*( this->MomentVector().size()/pow(2.0,30.0)+2.0*this->SystemSize()/pow(2.0,30.0) );
}
	
