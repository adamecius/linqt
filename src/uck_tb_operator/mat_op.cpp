
#include "uck_tb_operator.hpp"

void 
UnitCellK_TBOp::Rescale(const qt::real a, const qt::real b)
{
	for( qt::index i=0; i < val.size(); i++)
		val[i]=a*val[i];
	qt::real k[3]={0,0,0};
	Addk_dependence( k );
}

void UnitCellK_TBOp::Multiply(	const qt::real a,const qt::complex* X, 
				const qt::real b, qt::complex* Y) const 
{
	for( qt::index i=0; i < Dim(); i++)
	{
		qt::complex t= 0.0;
		for( qt::index j=0; j < Dim() ; j++)
		{
			const qt::complex O=Op(i,j);
			const qt::complex x=X[j];
			t.real(t.real() + x.real()*O.real() - x.imag()*O.imag()) ;
			t.imag(t.imag() + x.real()*O.imag() + x.imag()*O.real()) ;
		}
		qt::complex y=Y[i];
		y.real( b*y.real() + a*t.real());
		y.imag( b*y.imag() + a*t.imag());
		Y[i]= y;
	}	
}


void UnitCellK_TBOp::Print(	) const 
{
	if( Dim() > 10)
	for( qt::index i=0; i < Dim(); i++)
	{
		std::cout<<" Row = "<<i<<" : "<<std::endl;
		std::cout<<" Column\t Real \t Imag "<<std::endl;
		for( qt::index j=0; j < Dim(); j++)
		{
			std::cout<<j<<"\t "<<Op(i,j).real()<<"\t "<<Op(i,j).imag()<<std::endl;
		}
	}
	else
	for( qt::index i=0; i < Dim(); i++)
	{
		for( qt::index j=0; j < Dim(); j++)
			std::cout<<"("<<Op(i,j).real()<<","<<Op(i,j).imag()<<")\t ";
		std::cout<<std::endl;
	}
		
}
