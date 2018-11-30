#ifndef KPM_HPP
#define KPM_HPP

#include "types_definitions.hpp"
#include "kpm_config.hpp" //include parser and kpm::config structure
#include "sparse_block_matrix.hpp"


namespace kpm
{
	
	void cheb_evolve(	const qt::index n, 
						const sparse::BlockMatrix<qt::complex>& H,
						qt::complex* &jn0,
						qt::complex* &jn1
						);

	bool 
	GetOptionsFromCFG(std::ifstream& configFile, config& opt);

	void 
	PrintOptions(const config opt);


}



#endif
