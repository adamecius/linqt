
#include "uck_tb_operator.hpp"

bool UnitCellK_TBOp::readOperator( std::string inputName)
{
	std::ifstream inputFile;
	inputFile.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
	try  {	inputFile.open(inputName.c_str());	}
	catch (std::ifstream::failure e)
	{
		inputFile.close();
		std::cerr 	<<"Exception:"<<e.what()
					<<"\n while opening k-op  info sfile: "
					<< inputName<<"."<<std::endl;
			return false;
	}
	qt::integer 
		_OpDim,_numEntries; 
	qt::real 
		_latConst=0;	
	std::vector< std::vector < qt::real > >
		_lat(3); 
	for(int i=0;i<3; i++ )
		_lat[i] = std::vector < qt::real >(3);



	inputFile>>_OpDim>>numEntries;
	inputFile>>_latConst;
	inputFile>>_lat[0][0]>>_lat[0][1]>>_lat[0][2];
	inputFile>>_lat[1][0]>>_lat[1][1]>>_lat[1][2];
	inputFile>>_lat[2][0]>>_lat[2][1]>>_lat[2][2];


	SetDim(_OpDim);
	SetLatticeInfo(_latConst,_lat);

	Op = TBMat( Dim() , Dim() );
	val = std::vector< qt::complex > (numEntries);
	row = std::vector< qt::index > (numEntries);
	col = std::vector< qt::index > (numEntries);
	ishift = std::vector< std::vector<qt::real > > (numEntries);

	for(qt::index n=0; n<numEntries ; n++)
	{
		qt::real Reval,Imval;
		ishift[n] = std::vector<qt::real>(3);
		inputFile>>row[n]>>col[n]>>Reval>>Imval>>ishift[n][0]>>ishift[n][1]>>ishift[n][2];
		val[n] = qt::complex (Reval, Imval);
	}
	
	pos=ishift;
	for(qt::index n=0; n<numEntries ; n++)
	for(qt::index i=0; i<3 ; i++)
	{
		pos[n][i] = 0;
		for(qt::index m=0; m<3 ; m++)
			pos[n][i] += ishift[n][m]*lat[m][i];
	}

	qt::real k[3]={0,0,0};
	Addk_dependence( k );
};
