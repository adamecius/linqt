#ifndef HOPPING_KIND 
#define HOPPING_KIND

#include <random>
#include <complex>
#include <string>

using namespace std;

namespace hopping_kind_random
{
	static random_device rd;
};
//For defining dynamically the hopping kind, we will use the
//builder programming paradim. As it will allow for creating different 
//hopping instance, based on the users inputs.
class hopping_kind
{
	public:
	typedef complex<double> value_t;
	typedef default_random_engine random_generator;
	

	hopping_kind(value_t _hop): constant_hop(true), hop(_hop) { }; 

	hopping_kind(): 
	constant_hop(false),
	hop(0.0),
	biased_coin(0.0),
	rand_re(0.0,0.0),
	rand_im(0.0,0.0),
	gen(hopping_kind_random::rd())
	{
	}; 

	inline
	value_t operator()()  
	{ 
		if( constant_hop )
			return hop;

		if( biased_coin(gen) )
			return value_t(rand_re(gen),rand_im(gen));

		return  value_t(0,0);
	};
	
	inline		
	hopping_kind Concentration(double concentration)  
	{
		biased_coin = bernoulli_distribution(concentration);
		return *this;
	};

	inline		
	hopping_kind DistributionType(string& distribution,const value_t& min_val,const value_t& max_val ) 
	{
		this->min_r = min_val;
		this->max_r = max_val;
		rand_re		= uniform_real_distribution<>(min_val.real(),max_val.real());
		rand_im		= uniform_real_distribution<>(min_val.imag(),max_val.imag());
		return *this;
	};
	
	inline		
	void Rescale( const value_t& scalfact ) 
	{
		if( constant_hop ) 
			hop*= scalfact;
		else
			DistributionType(this->dist,scalfact*this->min_r,scalfact*this->max_r );
	};
	
	inline		
	bool is_constant() const 
	{
		return constant_hop;
	}

	private:
	const bool constant_hop;
	value_t hop, min_r, max_r;
	bernoulli_distribution biased_coin;
	uniform_real_distribution<>  rand_re,rand_im;
	random_generator gen;
	string dist;

};





#endif
