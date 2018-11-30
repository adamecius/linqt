
#include "uck_tb_operator.hpp"

void UnitCellK_TBOp::Addk_dependence( const qt::real* k )	
{
	qt::complex x=0;
	for( qt::index i=0; i < Dim(); i++)
	for( qt::index j=0; j < Dim(); j++)
		Op(i,j) = 0.0;

	for( qt::index n=0; n < numEntries ; n++ )
	{
		const qt::real dotv= pos[n][0]*k[0] + pos[n][1]*k[1] + pos[n][2]*k[2];
		qt::complex expv= qt::complex( cos(dotv) , -sin(dotv) );
		qt::complex v= val[n];
		qt::complex O= qt::complex( v.real()*expv.real() - v.imag()*expv.imag() , v.real()*expv.imag() + v.imag()*expv.real() );
		Op(row[n],col[n])+=O;
	}
}
