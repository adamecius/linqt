
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

	UnitCellK_TBOp(const qt::integer _OpDim=1):OpDim_(_OpDim){};	

	UnitCellK_TBOp(std::string inputName)
	{
		readOperator(inputName);
	};
	// PUBLIC COMMON FUNCTIONS
	public: 

	//
	void AddEntry( 	const qt::integer r, //row
					const qt::integer c, //column
					const qt::complex v, //value
					const qt::real is0,//shift in i0 
					const qt::real is1,//shift in i1 
					const qt::real is2);//shift in i2

	// This function add the k-depence
	void Addk_dependence( const qt::real* k );


	void writeOperator( std::string outputname );
	
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

	inline void
	SetLatticeInfo(const qt::real _latConst, const std::vector< std::vector < qt::real > >_lat)
	{
		lat=_lat; rec=_lat; //initialize for the same size
		latconst = _latConst ;
		for(qt::index m=0; m<3 ; m++)
		for(qt::index n=0; n<3 ; n++)
			lat[m][n] = _latConst *lat[m][n];
		
		latVol =-lat[0][2]*lat[1][1]*lat[2][0] + 
				 lat[0][1]*lat[1][2]*lat[2][0] +
				 lat[0][2]*lat[1][0]*lat[2][1] - 
				 lat[0][0]*lat[1][2]*lat[2][1] - 
				 lat[0][1]*lat[1][0]*lat[2][2] +
				 lat[0][0]*lat[1][1]*lat[2][2];

	rec =_lat ; //Initialize with the same size
	rec[0][0]=-lat[1][2]*lat[2][1] + lat[1][1]*lat[2][2]; 
	rec[0][1]= lat[1][2]*lat[2][0] - lat[1][0]*lat[2][2]; 
	rec[0][2]=-lat[1][1]*lat[2][0] + lat[1][0]*lat[2][1];

	rec[1][0]= lat[0][2]*lat[2][1] - lat[0][1]*lat[2][2]; 
	rec[1][1]=-lat[0][2]*lat[2][0] + lat[0][0]*lat[2][2]; 
	rec[1][2]= lat[0][1]*lat[2][0] - lat[0][0]*lat[2][1];

	rec[2][0]=-lat[0][2]*lat[1][1] + lat[0][1]*lat[1][2]; 
	rec[2][1]= lat[0][2]*lat[1][0] - lat[0][0]*lat[1][2]; 
	rec[2][2]=-lat[0][1]*lat[1][0] + lat[0][0]*lat[1][1];

	for(qt::index m=0; m<3 ; m++)
	for(qt::index n=0; n<3 ; n++)
		rec[m][n]=rec[m][n]*2.0*M_PI/latVol;
		  };
	
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
	std::vector< std::vector<qt::real> > ishift;
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
