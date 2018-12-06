
#include "uck_tb_operator.hpp"

void 
UnitCellK_TBOp::writeOperator( std::string output_name)
{
	output_name+=".KOP";
	std::ofstream output( output_name.c_str() );

	
	qt::integer numEntry=0.0;
	for(qt::integer n=0; n< val.size(); n++)
	if( std::norm(val[n])!=0 )
			numEntry++;
	
	output<<Dim()<<" "<<numEntry<<" "<<std::endl;
	for(qt::integer n=0; n< val.size(); n++)
		if( std::norm(val[n])!=0 )
			output	<<row[n]<<" "<<col[n]<<" "
					<<val[n].real()<<" "<<val[n].imag()<<" "
					<<ishift[n][0]<<" "<<ishift[n][1]<<" "<<ishift[n][2]<<std::endl;		
		output.close();
};
