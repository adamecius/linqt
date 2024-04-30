#include "chebyshev_moments.hpp"



void chebyshev::Moments::SetInitVectors( SparseMatrixType &OP ,const Moments::vector_t& T0 )
{
	const auto dim = this->SystemSize();
	assert( OP.rank() == this->SystemSize() && T0.size() == this->SystemSize() );

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	linalg::copy ( T0, this->Chebyshev1() );
	OP.Multiply( this->Chebyshev1(), this->Chebyshev0() );
	this->Hamiltonian().Multiply( this->Chebyshev0(), this->Chebyshev1() );

	return ;
};



	//light functions
int chebyshev::Moments::JacksonKernelMomCutOff( const double broad )
{
	assert( broad >0 );
	const double eta   =  2.0*broad/1000/this->BandWidth();
	return ceil(M_PI/eta);
};
	
//light functions
double chebyshev::Moments::JacksonKernel(const double m,  const double Mom )
{
	const double
	phi_J = M_PI/(double)(Mom+1.0);
	return ( (Mom-m+1)*cos( phi_J*m )+ sin(phi_J*m)*cos(phi_J)/sin(phi_J) )*phi_J/M_PI;
};




void chebyshev::Moments::SetInitVectors( const Moments::vector_t& T0 )
{
	assert( T0.size() == this->SystemSize() );
	const auto dim = this->SystemSize();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	linalg::copy ( T0, this->Chebyshev0() );
	this->Hamiltonian().Multiply( this->Chebyshev0(), this->Chebyshev1() );
};



void chebyshev::Moments1D_nonOrth::SetInitVectors_nonOrthogonal( Moments::vector_t& T0 )
{
  

	assert( T0.size() == this->SystemSize() );
	const auto dim = this->SystemSize();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	
        vector_t tmp_(dim,0.0);



        linalg::orthogonalize(*S_, this->Chebyshev0(), T0);
	
	
	this->Hamiltonian().Multiply( this->Chebyshev0(), tmp_ );
        linalg::orthogonalize(*S_, this->Chebyshev1(), tmp_);
	
};


void chebyshev::Moments1D_nonOrth::SetInitVectors_nonOrthogonal( SparseMatrixType &OP , Moments::vector_t& T0 )
{
	const auto dim = this->SystemSize();
	assert( OP.rank() == this->SystemSize() && T0.size() == this->SystemSize() );

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()


	vector_t tmp_(dim,0.0);


	linalg::copy ( T0, this->Chebyshev1() );
	OP.Multiply( this->Chebyshev1(), this->Chebyshev0() );
	this->Hamiltonian().Multiply( this->Chebyshev0(), this->Chebyshev1() );



	

	linalg::orthogonalize(*S_, tmp_, T0);
	linalg::copy ( tmp_, this->Chebyshev1() );
        OP.Multiply( this->Chebyshev1(), this->Chebyshev0() );

	
	this->Hamiltonian().Multiply( this->Chebyshev0(), tmp_ );
        linalg::orthogonalize(*S_, this->Chebyshev1(), tmp_ );


		
	return ;
};


int chebyshev::Moments1D_nonOrth::Iterate_nonOrthogonal( )
{

        std::size_t dim = Chebyshev0().size();
        vector_t tmp_( dim, 0.0 ),
	         tmp_2_( dim, 0.0);


	
	this->Hamiltonian().Multiply( this->Chebyshev1(),tmp_);
	
       	linalg::orthogonalize(*S_, tmp_2_, tmp_);
	//linalg::copy(tmp_, tmp_2_);
	
	#pragma omp parallel for
	for(std::size_t i=0; i< dim;i++){
	  Chebyshev0()[i] = 2.0 * tmp_2_[i] - Chebyshev0()[i];
	}
	
	this->Chebyshev0().swap(this->Chebyshev1());

	

	return 0;
};



int chebyshev::Moments::Iterate( )
{
	this->Hamiltonian().Multiply(2.0,this->Chebyshev1(),-1.0,this->Chebyshev0());
	this->Chebyshev0().swap(this->Chebyshev1());
	return 0;
};
