
#include "uck_tb_operator.hpp"

void UnitCellK_TBOp::readLatticeInfo( std::string inputName)
{
	lat = std::vector< std::vector< qt::real> >(3);
	rec = std::vector< std::vector< qt::real> >(3);
	for(qt::index i=0;i<3;i++)
	{
		lat[i] = std::vector< qt::real>(3);
		rec[i] = std::vector< qt::real>(3);
	}

	std::ifstream inputFile;
	inputFile.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
	try  {	inputFile.open(inputName.c_str());	}
	catch (std::ifstream::failure e)
	{
		inputFile.close();
		std::cerr 	<<"Exception:"<<e.what()
					<<"\n while opening lattice info sfile: "
					<< inputName<<"."<<std::endl;
		std::exit(-1);
	}
		
	inputFile>>latconst;
		
	for(qt::index m=0; m<3 ; m++)
	for(qt::index n=0; n<3 ; n++)
	{
		inputFile>>lat[m][n];
		lat[m][n]=lat[m][n]*latconst;
	}
	inputFile.close();
		
	latVol =-lat[0][2]*lat[1][1]*lat[2][0] + 
			 lat[0][1]*lat[1][2]*lat[2][0] + 
			 lat[0][2]*lat[1][0]*lat[2][1] - 
			 lat[0][0]*lat[1][2]*lat[2][1] - 
			 lat[0][1]*lat[1][0]*lat[2][2] +
			 lat[0][0]*lat[1][1]*lat[2][2];


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
}
