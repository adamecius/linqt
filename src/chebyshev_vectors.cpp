#include "chebyshev_moments.hpp"


void chebyshev::Vectors::SetInitVectors( SparseMatrixType &NHAM,const Moments::vector_t& T0 )
{
	assert( NHAM.rank() == this->SystemSize()&& T0.size() == this->SystemSize() );
	
	if( ChebV0.size()!= T0.size() )
		ChebV0 = Moments::vector_t(T0.size(),Moments::value_t(0)); 

	if( ChebV1.size()!= T0.size() )
		ChebV1 = Moments::vector_t(T0.size(),Moments::value_t(0)); 
	
	linalg::copy( T0, ChebV0 );
	NHAM.Multiply( ChebV0, ChebV1 );
};


void chebyshev::Vectors::SetInitVectors( SparseMatrixType &NHAM, SparseMatrixType &OP ,const Moments::vector_t& T0 )
{
	assert( OP.rank() == NHAM.rank() && NHAM.rank() == this->SystemSize()&& T0.size() == this->SystemSize() );

	if( ChebV0.size()!= T0.size() )
		ChebV0 = Moments::vector_t(T0.size(),Moments::value_t(0)); 

	if( ChebV1.size()!= T0.size() )
		ChebV1 = Moments::vector_t(T0.size(),Moments::value_t(0)); 

	linalg::copy ( T0, ChebV1 );
	OP.Multiply  ( ChebV1, ChebV0 );
	NHAM.Multiply( ChebV0, ChebV1 );
};


int chebyshev::Vectors::IterateAll( SparseMatrixType &NHAM )
{
	const size_t  dim = NHAM.rank();
	assert( dim == this->SystemSize()  );
	
	linalg::copy( this->ChebV0 ,this->Vector(0) );
	for(int m=1; m < this->HighestMomentNumber(); m++ )
	{
		linalg::copy( ChebV1 , this->Vector(m) );
		NHAM.Multiply(2.0,ChebV1,-1.0,ChebV0);
		ChebV0.swap(ChebV1);
	}
	return 0;
};

int chebyshev::Vectors::EvolveAll( SparseMatrixType &NHAM, const double DeltaT, const double Omega0)
{
	const auto dim = NHAM.rank();
	const auto numMom = this->HighestMomentNumber();

	if( this->ChebV0.size()!= dim )
		this->ChebV0 = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->ChebV1.size()!= dim )
		this->ChebV1 = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in ChebV0 and ChebV1

	const auto I = Moments::value_t(0, 1);
	const double x = Omega0*DeltaT;
	for(size_t m=0; m < numMom; m++ )
	{
		auto& myVec = this->Vector(m);
		
		int n = 0;
		double Jn = besselJ(n,x);
		linalg::copy(myVec , ChebV0);
		linalg::scal(0, myVec); //Set to zero
		linalg::axpy( Jn , ChebV0, myVec);

		double Jn1 = besselJ(n+1,x);	
		NHAM.Multiply(ChebV0, ChebV1);
		linalg::axpy(-2 * I * Jn1, ChebV1, myVec);
		
		auto nIp =-I;
		while( 0.5*(std::abs(Jn)+std::abs(Jn1) ) > 1e-15)
		{
			nIp*=-I ;
			Jn  = Jn1;
			Jn1 = besselJ(n, x);
			NHAM.Multiply(2.0, ChebV1, -1.0, ChebV0);
			linalg::axpy(2 * nIp *Jn1, ChebV0, myVec);
			ChebV0.swap(ChebV1);
			n++;
		}
	}
  return 0;
};

int chebyshev::Vectors::Multiply( SparseMatrixType &OP )
{
	assert( OP.rank() == this->SystemSize() );
	const auto dim = OP.rank();
	const auto numMom = this->HighestMomentNumber();
	vectorList_t*  pChebmu = &(this->Chebmu);
	const int nthreads = mkl_get_max_threads();
	auto  pOPV 	  = &this->OPV;
		
	if( this->OPV.size()!= dim )
		this->OPV = Moments::vector_t ( dim );
	
	for(size_t m=0; m < numMom; m++ )
	{
		linalg::copy( pChebmu->ListElem(m), *pOPV ); 
		OP.Multiply( *pOPV, pChebmu->ListElem(m) );
	}

	return 0;
};

double chebyshev::Vectors::MemoryConsumptionInGB()
{
	return sizeof(value_t)*( this->Size()/pow(2.0,30.0)+2.0*this->SystemSize()/pow(2.0,30.0) );
}
	
