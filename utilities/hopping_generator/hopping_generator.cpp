
#include "lattice_geometry.h"
#include <fstream>
#include <string>
#include <complex>
#include <iostream>
int main(int argc, char *argv[])
{
	if( argc != 2)
	{
		std::cout<<"submit a name for the output"<<std::endl;
		return 0;
	}
	std::string hopOutputFname(argv[1]);
	std::ofstream hopFile((hopOutputFname+".hop").c_str() );

	const int
	i0=0, i1=0, i2=0 ;

//Ws2
//	Complex E0[2]   ={ -1.31/1000. , 1.31/1000. };
//	Complex t1      = -2.657;
//	Complex t2      =0. ;
//	Complex VI[2]   ={-1.02/1000. , 1.21/1000.	 };
//	Complex VR      = 0.36/1000. ;
//	Complex VPIA[2] = { -0.98/1000., -3.81};

//Wse2
//	Complex E0[2]   ={ -0.54/1000. , 0.54/1000. };
//	Complex t1      = -2.507;
//	Complex t2      =0. ;
//	Complex VI[2]   ={-1.22/1000. , 1.16/1000.	 };
//	Complex VR      = 0.56/1000. ;
//	Complex VPIA[2] = { -2.69/1000., -2.54/1000.};

//PIA
	Complex E0[2]   ={ 0. , 0. };
	Complex t1      = -2.8;
	Complex t2      =0. ;
	Complex VI[2]   ={ 0. , 0.	 };
	Complex VR      = 0. ;
	Complex VPIA[2] = { 0.2, 0.2 };




	for(int io=0; io < 2 ; io++ )
		for(int  is=0; is < 2 ; is++ )
		{
			const Real	
			ss= 1 - 2*is,
			oo= 1 - 2*io ;

			const int
			jo=1-io,
			js=1-is,
			dI=1-2*io;

			int J0nn[3]={i0,i0+dI,i0   };
			int J1nn[3]={i1,i1   ,i1+dI};

			int J0nnn[6]={i0-1,i0+1,i0+0,i0+0,i0-1,i0+1};
			int J1nnn[6]={i1+0,i1+0,i1-1,i1+1,i1+1,i1-1};

//----------------------------ONSITE----------------------------/
			if( std::norm(E0[io]) != 0)
			hopFile<<io<<" "<<is<<" "
				   <<i0<<" "<<i1<<" "<<i2<<" "
				   <<io<<" "<<is<<" "
				   <<E0[io].real()<<" "<<E0[io].imag()<<std::endl;

//----------------------------NEAREST NEIGHBORS/----------------------------/
			if( std::norm(t1) != 0)
			for(int k=0;k<3;++k)
				hopFile<<io<<" "<<is<<" "
					   <<J0nn[k]<<" "<<J1nn[k]<<" "<<i2<<" "
					   <<jo<<" "<<is<<" "
					   <<t1.real()<<" "<<t1.imag()<<std::endl;

//----------------------------NEXT NEAREST NEIGHBORS/----------------------------/
			if( std::norm(t2)!=0 || std::norm(VI[io])!=0)
			{
				Complex val[6]=
				{
						(t2- I*ss*oo*VI[io]/3./sqrt(3) ),
						(t2+ I*ss*oo*VI[io]/3./sqrt(3)),
						(t2+ I*ss*oo*VI[io]/3./sqrt(3)),
						(t2- I*ss*oo*VI[io]/3./sqrt(3)),
						(t2+ I*ss*oo*VI[io]/3./sqrt(3)),
						(t2- I*ss*oo*VI[io]/3./sqrt(3))
				};
				for(int k=0;k<6;++k)
					hopFile<<io<<" "<<is<<" "
						   <<J0nnn[k]<<" "<<J1nnn[k]<<" "<<i2<<" "
						   <<io<<" "<<is<<" "
						   <<val[k].real()<<" "<<val[k].imag()<<std::endl;
			}
//---------------------------- NEAREST NEIGHBORS SPIN FLIP/----------------------------/
			if( std::norm(VR)!=0 )
			{
				const Real 
				ri[3]={ i0*A[0][0] + i1*A[1][0] + io*Delta[0] ,
						i0*A[0][1] + i1*A[1][1] + io*Delta[1] , 
						0 };

				for(int k=0;k<3;++k)
				{

				const Real 
				rj[3]={ J0nn[k]*A[0][0] + J1nn[k]*A[1][0] + jo*Delta[0] ,
					J0nn[k]*A[0][1] + J1nn[k]*A[1][1] + jo*Delta[1] ,
						0 
						};
						
				const Real 
				del_r[3]={ ri[0]-rj[0] ,ri[1]-rj[1] ,ri[2]-rj[2]  };

				const Complex PhiRSO= CrossProductDotZ(del_r, is,js );
			//	std::cout<<PhiRSO.real()<<" "<<PhiRSO.imag()<<std::endl;

					Complex val=2.*I*VR*PhiRSO/3.;
					hopFile<<io<<" "<<is<<" "
						   <<J0nn[k]<<" "<<J1nn[k]<<" "<<i2<<" "
						   <<jo<<" "<<js<<" "
						   <<val.real()<<" "<<val.imag()<<std::endl;
				}
			}

			//------ PIA INTERACTION------
			if( std::norm(VPIA[io]) !=0)
			{

				const Real 
				ri[3]={ i0*A[0][0] + i1*A[1][0] + io*Delta[0] ,
						i0*A[0][1] + i1*A[1][1] + io*Delta[1] , 
						0 };

				for(int k=0;k<6;++k)
				{

					const Real 
					rj[3]={ J0nnn[k]*A[0][0] + J1nnn[k]*A[1][0] + io*Delta[0] ,
						J0nnn[k]*A[0][1] + J1nnn[k]*A[1][1] + io*Delta[1] ,
						0 
							};
						
					const Real 
					del_r[3]={ ri[0]-rj[0] ,ri[1]-rj[1] ,ri[2]-rj[2]  };

					const Complex PhiPIA= CrossProductDotZ(del_r, is,js );

					const Complex 
					val= 2.*I*VPIA[io]*PhiPIA/3.;
					hopFile<<io<<" "<<is<<" "
						   <<J0nnn[k]<<" "<<J1nnn[k]<<" "<<i2<<" "
						   <<io<<" "<<js<<" "
						   <<val.real()<<" "<<val.imag()<<std::endl;
				}
			}
			
		}
	hopFile.close();
return 0;}

