#ifndef KPM_HPP
#define KPM_HPP

#include "chebyshev_set.hpp"
#include "kpm_linalg.hpp"

#include "kpm_utilities.hpp"
#include "types_definitions.hpp"
#include <iostream>
#include "random.hpp"
#include "lattice.hpp"


bool IsNumber(const std::string& s);



///Maybe we need to move it to a utility file
NumCal::integer
GetMomentsPerNode( const NumCal::integer numMoms);

/**In this class we will include everything needed to perform
 * a calculation using the Kernel Polynomial Method*/

namespace NumCal{

/** \addtogroup kpm
 *  @{
 */

/**\class Kpm
 *
 * \brief Set of approximation based on KPM
 *
 * This class is aimed to be used as a container of the different sets
 * of approximations done within the Kernel Polynomial method approximation
 *
 *
 * \note This class is still in development process, use it at your own risk!
 *
 * \version 1.0
 *
 */
class Kpm
{
	public:

	/********************************CONSTRUCTORS*****************************/
	///The Main constructor
	/*!
	 * This constructor will receive the dimension
	 * of the working Hilber space and will try to
	 * open the file _config_filename. If succeed
	 * will read the following parameters:\n
	 * Truncation Order.\n
	 * Minimal Bound\n
	 * Maximal Bound\n
	 * CutOff\n
	 */
	Kpm(const NumCal::integer _dim,const  std::string _config_filename):
	vecDim_(_dim),
	chebSet( ChebyshevSet("", _dim) ),
	linalg (_dim, chebSet.MemSep() )
	{
		std::ifstream config_file(_config_filename.c_str());
		///Look for the line kpm_infor
		std::string line;
		bool found=false;
		for (int line_num = 0; std::getline(config_file, line); ++line_num)
			if(line=="kpm_info")
			{
				found=true;
				break;
			}
		if(found)
		{
			NumCal::integer trunc_order;
			NumCal::real min_bound, max_bound,cutoff;
			config_file>>trunc_order>>min_bound>>max_bound>>cutoff;

			seed_=234324;


			SetTruncOrder(trunc_order);
			SetBounds(min_bound,max_bound,cutoff);
			config_file.close();
		}else
		{
			config_file.close();
			std::cerr<<"The config file does not posses a kpm_info section. The simulation cannot proceed"<<std::endl;
			std::exit(-1);
		}
		distribution= NumCal::random::uniform_real_dist( -sqrt(3) ,  sqrt(3) ) ;
	}
/***********************Convertion functions**************************/
	NumCal::real EnormToEn(const NumCal::real Enorm) const;

	NumCal::real EnToEnorm(const NumCal::real En) const;


/**************************SETTERS*************************************/
	void SetTruncOrder(const NumCal::size_t _trunc_order);

	void SetBounds(	const NumCal::real _min_bound,
					const NumCal::real _max_bound,
					const NumCal::real _cutoff);

/**************************GETTERS*************************************/
	NumCal::size_t TruncOrder() const
	{
		return trunc_order_;
	};

	NumCal::real MinBound() const
	{
		return min_bound_;
	};

	NumCal::real MaxBound() const
	{
		return max_bound_;
	};

	NumCal::real CutOff() const
	{
		return cutoff_;
	};

	NumCal::real EnergyFactor() const
	{
		return scale_factor_;
	};

	NumCal::integer VecDim() const
	{
		return vecDim_;
	}

	void SetRandomSeed(unsigned int _seed)
	{
		seed_=_seed;
		generator.seed( seed_ );
	}

	unsigned int RandomSeed()
	{
		return seed_;
	}

	void ConductivityMoments(Lattice _lattice, std::string _condMomXXFile,std::string _condMomXYFile);

	void SpinConductivityMoments(Lattice _lattice, std::string _condMomXXFile,std::string _condMomXYFile);

	void ValleyConductivityMoments(const NumCal::integer vindex,const NumCal::real RKscal,  Lattice _lattice, std::string
_condMomXXFile,std::string 
_condMomXYFile);

	void DOSMoments(Lattice _lattice, std::string _dosMomFile);

	/**************************PRIVATE VARIABLES*************************************/

	private:
		NumCal::integer vecDim_, trunc_order_;
		NumCal::real	min_bound_,max_bound_,cutoff_,scale_factor_;
		NumCal::random::uniform_real_dist distribution;
		NumCal::random::generator generator;
		unsigned int seed_;

		ChebyshevSet chebSet;
		KPMLinalg linalg;
};

/** @}*/

}

#endif
