#include "chebyshev_moments.hpp"


int chebyshev::Vectors::IterateAll( )
{	
	linalg::copy( this->Chebyshev0() ,this->Vector(0) );
	for(int m=1; m < this->HighestMomentNumber(); m++ )
	{
		linalg::copy( Chebyshev1() , this->Vector(m) );
		this->Hamiltonian().Multiply(2.0,Chebyshev1(),-1.0,Chebyshev0());
		Chebyshev0().swap(Chebyshev1());
	}
	return 0;
};

int chebyshev::Vectors::EvolveAll(const double DeltaT, const double Omega0)
{
	const auto dim = this->SystemSize();
	const auto numMom = this->HighestMomentNumber();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	const auto I = Moments::value_t(0, 1);
	const double x = Omega0*DeltaT;
	for(size_t m=0; m < numMom; m++ )
	{
		auto& myVec = this->Vector(m);
		
		int n = 0;
		double Jn = besselJ(n,x);
		linalg::copy(myVec , Chebyshev0());
		linalg::scal(0, myVec); //Set to zero
		linalg::axpy( Jn , Chebyshev0(), myVec);

		double Jn1 = besselJ(n+1,x);	
		this->Hamiltonian().Multiply(Chebyshev0(), Chebyshev1());
		linalg::axpy(-2 * I * Jn1, Chebyshev1(), myVec);
		
		auto nIp =-I;
		while( 0.5*(std::abs(Jn)+std::abs(Jn1) ) > 1e-15)
		{
			nIp*=-I ;
			Jn  = Jn1;
			Jn1 = besselJ(n, x);
			this->Hamiltonian().Multiply(2.0, Chebyshev1(), -1.0, Chebyshev0());
			linalg::axpy(2 * nIp *Jn1, Chebyshev0(), myVec);
			Chebyshev0().swap(Chebyshev1());
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
	return SizeInGB()+2.0*( (double)this->SystemSize() )*pow(2.0,-30.0) ;
}
	
