
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
	qt::integer OpDim,_numEntries; 
	inputFile>>OpDim>>numEntries;
	SetDim(OpDim);

	Op = TBMat( Dim() , Dim() );
	val = std::vector< qt::complex > (numEntries);
	row = std::vector< qt::index > (numEntries);
	col = std::vector< qt::index > (numEntries);
	pos = std::vector< std::vector<qt::real > > (numEntries);
	
	for(qt::index n=0; n<numEntries ; n++)
		pos[n] = std::vector<qt::real>(3);
		
	for(qt::index n=0; n<numEntries ; n++)
	{
		qt::real Reval,Imval;
		inputFile>>row[n]>>col[n]>>Reval>>Imval>>pos[n][0]>>pos[n][1]>>pos[n][2];
		val[n] = qt::complex (Reval, Imval);
	}
	
	qt::real k[3]={0,0,0};
	Addk_dependence( k );
};
