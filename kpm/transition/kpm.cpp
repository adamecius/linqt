

#include "kpm.hpp"


my::real
Kpm::EnormToEn(const my::real Enorm) const
{
	return EnergyFactor()*Enorm+0.5*(MaxBound()+MinBound());
};

my::real
Kpm::EnToEnorm(const my::real En) const
{
	return (En-0.5*(MaxBound()+MinBound()) )/EnergyFactor();
};

/**********************KERNEL FUNCTIONS********************************/
my::real
Kpm::jackson_kernel(const my::real m)
{
const my::real M=TruncOrder();
const my::real THETA=M_PI/(M+1);
return
((M-m+1)*cos(m*THETA)+sin(m*THETA)*cos(THETA)/sin(THETA))/(M+1);
}
/**************************SETTERS*************************************/
void
Kpm::SetTruncOrder(const my::size_t _trunc_order)
{
	if( _trunc_order <= 0 )
	{
		std::cerr<<" The trucation order :"<<_trunc_order;
		std::cerr<<" should be higher than zero"<<std::endl;
		std::exit(-1);
	}

	trunc_order_=_trunc_order;
};

void
Kpm::SetBounds(	const my::real _min_bound,
				const my::real _max_bound,
				const my::real _cutoff)
{
	if( _max_bound <= _min_bound)
	{
		std::cerr<<"The bounds:"<<_max_bound<<" "<<_min_bound;
		std::cerr<<" are incorrect"<<std::endl;
		std::exit(-1);
	}
	if( _cutoff >= 1 )
	{
		std::cerr<<"The cutoff:"<<_cutoff;
		std::cerr<<" is incorrect"<<std::endl;
		std::exit(-1);
	}
	min_bound_=_min_bound;
	max_bound_=_max_bound;
	cutoff_=_cutoff;
	scale_factor_= (MaxBound()-MinBound())/( 2.*CutOff() );
};

void
Kpm::SetEnergyRange( const my::real _Emin,
					 const my::real _Emax,
					 const my::size_t _NE)
{
	if( _Emin < MinBound() )
	{
		std::cerr<<"The minimum energy range:"<<_Emin<<" is smaller";
		std::cerr<<" than the minimal bound"<<MinBound()<<std::endl;
		std::exit(-1);
	}
	if( _Emax > MaxBound() )
	{
		std::cerr<<"The maximum energy range:"<<_Emax<<" is smaller";
		std::cerr<<" than the maximum bound"<<MaxBound()<<std::endl;
		std::exit(-1);
	}
	if( _NE <= 0 )
	{
		std::cerr<<" The number of energy points:"<<_NE;
		std::cerr<<" should be higher than zero"<<std::endl;
		std::exit(-1);
	}

	Emin_=_Emin;
	Emax_=_Emax;
	NE_=_NE;
	dE_= (Emax_-Emin_)/_NE;
};




