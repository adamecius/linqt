
#ifndef UCK_TB_OP_HPP
#define UCK_TB_OP_HPP

#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include "types_definitions.hpp"
#include "dense_matrix.hpp" 


class UnitCellK_TBOp
{
	typedef qt::dense_matrix<qt::complex> TBMat;

	// THE CONSTRUCTORS
	public:

	UnitCellK_TBOp(void){};	

	UnitCellK_TBOp(std::string inputName)
	{
		readOperator(inputName);
	};
	// PUBLIC COMMON FUNCTIONS
	public: 

	// This function add the k-depence
	void Addk_dependence( const qt::real* k );

	// This function read the operator from a file name inputName
	bool readOperator( std::string inputName);

	// This function read the lattice from a file name inputName
	void readLatticeInfo( std::string inputName);	


	//This function perform a multiplication of the matrix
	void Print() const;


	// This function rescale the matrix
	void 
	Rescale(const qt::real a, const qt::real b=0.0);

	//This function perform a multiplication of the matrix
	void Multiply(	const qt::real a,const qt::complex* X, 
					const qt::real b, qt::complex* Y) const;

	// PUBLIC GETTERS FUNCTIONS
	public: 

	inline void
	SetDim(const qt::dimension _OpDim){ OpDim_=_OpDim;  }
	
	inline qt::dimension
	Dim()  const { return OpDim_;  }
	
	inline TBMat
	Mat()  const { return Op;  }
	
	inline qt::complex 
	operator()( const qt::index i, qt::index j) const { return Op(i,j); }

	inline qt::real
	Lat( const qt::index i, const qt::index j) const { return lat[i][j]; }

	inline qt::real 
	LatConst() const { return latconst; }

	inline qt::real
	Rec( const qt::index i, const qt::index j) const { return rec[i][j]; }

	inline qt::index LatVol() const {return latVol;}

// Private variables
	private:
	qt::dimension OpDim_,numEntries;
	TBMat Op;
	std::vector< qt::complex > val;
	std::vector< qt::index > row, col;
	std::vector< std::vector<qt::real> > pos;
	std::vector< std::vector<qt::real> > lat;
	std::vector< std::vector<qt::real> > rec;
	double latVol, latconst;

};
/*
#include "uck_tb_operator/read_operator.cpp"
#include "uck_tb_operator/read_lattice.cpp"
#include "uck_tb_operator/add_k_dep.cpp"
#include "uck_tb_operator/mat_op.cpp"
*/
#endif
