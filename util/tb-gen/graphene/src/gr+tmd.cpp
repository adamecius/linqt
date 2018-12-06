
#include <fstream>
#include <string>
#include <iostream>
#include "types_definitions.hpp"
#include "parser.hpp"
#include "uck_tb_operator.hpp"


qt::complex
SIGMA[3][2][2]= 
{
  { { 0 , 1 },{ 1 , 0 } }	,
  { { 0 ,qt::complex(0,-1) },{ qt::complex(0,1) , 0 } }	,
  { { 1 , 0 },{ 0 ,-1 } }	
};

qt::complex 
CrossProductDotZ(const qt::real x[3], const qt::integer s0, const qt::integer s1  )
{
		qt::real norm= sqrt( x[0]*x[0] + x[1]*x[1]  ) ; //normalization vector
		return (SIGMA[0][s0][s1]*x[1] - SIGMA[1][s0][s1]*x[0])/norm;
};


int main(int argc, char *argv[])
{

	//Set the system label and try to open the configuration file
	// If fail to open it, abort the program
	std::string configFilename; 
	if ( argc!= 1) configFilename= argv[1];
	std::ifstream configFile( configFilename.c_str() );
	try{ 
		configFile.exceptions(configFile.failbit);
	}catch (const std::ios_base::failure& e){
		std::cerr 	<<"Exception:"<<e.what()
					<<"\n while opening config file: "
					<< configFilename <<std::endl;
		return 0;
	}
	

	std::string label, suffix, header;

	parser::GetValue(configFile,"Label" , label, parser::Optional);
	parser::GetValue(configFile,"Suffix", suffix,parser::Optional);


	qt::integer num_orb=0, max_spin=0;
	qt::real hop=0, StdPot;
	header="LATTICE_INFO";
	parser::GetValue(configFile,"NumberOfOrbitals", num_orb);
	parser::GetValue(configFile,"MaxSpin", max_spin);
	parser::GetValue(configFile,"HoppingEnergy", hop);
	parser::GetValue(configFile,"StaggeredPot", StdPot);
	const int orbsPerCell= num_orb*max_spin;


	std::vector< std::vector<qt::real> > Lat;
 	if(!parser::GetBlock(configFile,"LAT_VEC", Lat,3,3,parser::NotOptional  ))
 	{
		std::cerr 	<<"Problem Reading LAT_VEC block"<<std::endl;
		return 0;		
	}

	std::vector< std::vector<qt::real> > Opos;
 	if(!parser::GetBlock(configFile,"ORB_POS", Opos,num_orb,3,parser::NotOptional  ))
 	{
		std::cerr 	<<"Problem Reading ORB_POS block"<<std::endl;
		return 0;		
	}

			
	std::cout	<<"\nCreating the hamiltonian for a lattice with:"<<std::endl
				<<"with NumberOfOrbitals:\t" << num_orb<<std::endl
				<<"and maximum spin:\t"<<max_spin<<std::endl
				<<"The leading hopping energy is :\t"<<hop<<" eV "<<std::endl
				<<"The staggered potential is :\t"<<StdPot<<" eV "<<std::endl;

	std::cout<<"The system will be created using the following lattice vectors"<<std::endl;
	for(int i=0;i<3;i++)
	{
		std::cout<<"Lat_"<<i<<" = ";
		for(int j=0;j<3;j++)
				std::cout<<" "<<Lat[i][j];
		std::cout<<std::endl;
	}

	std::cout<<"The system has orbitals which are located at the following positions"<<std::endl;
	for(int i=0;i<num_orb;i++)
	{
		std::cout<<"OrbPos_"<<i<<" = ";
		for(int j=0;j<3;j++)
				std::cout<<" "<<Opos[i][j];
		std::cout<<std::endl;
	}


	qt::real RSO=0, ISOA=0, ISOB=0, PIAA=0, PIAB=0;
	header="SOC_INFO";
	if (parser::FindHeader(configFile, header )!= std::string::npos && max_spin==2 )
	{
		std::cout<<"Including spin-orbit exchange"<<std::endl;
		parser::GetValue(configFile,"Rashba", RSO);
		parser::GetValue(configFile,"IntrinsicA", ISOA);
		parser::GetValue(configFile,"IntrinsicB", ISOB);
		parser::GetValue(configFile,"PIA_A", PIAA);
		parser::GetValue(configFile,"PIA_B", PIAB);
		std::cout	<<" Rashba SOC = :\t"<<RSO<<" eV "<<std::endl
					<<" Rashba IntrinsicA = :\t"<<ISOA<<" eV "<<std::endl
					<<" Rashba IntrinsicB = :\t"<<ISOB<<" eV "<<std::endl
					<<" Rashba PIA_A = :\t"<<PIAA<<" eV "<<std::endl
					<<" Rashba PIA_B = :\t"<<PIAB<<" eV "<<std::endl;
	}


	header="EXCHANGE_INFO";
	qt::real Jex, Angex[2], Jtheta, Jphi,Jex_x,Jex_y,Jex_z;
	if (parser::FindHeader(configFile, header )!= std::string::npos && max_spin==2 )
	{
		std::cout<<"Including magnetic exchange"<<std::endl;
		parser::GetValue(configFile,"ExchangeCoupling", Jex,parser::NotOptional);
		parser::GetArray(configFile,"ExchangeAngles", &Angex[0],2,parser::NotOptional); 
		Jphi=Angex[0]*M_PI/180.; Jtheta=Angex[1]*M_PI/180.;
		Jex_x= Jex*cos(Jphi)*sin(Jtheta);
		Jex_y= Jex*sin(Jphi)*sin(Jtheta);
		Jex_z= Jex*cos(Jtheta);

		std::cout<<"with exchange vector (meV):" <<Jex<<" "<<std::endl;
		std::cout<<"with exchange vector (meV): ( " 
						<<Jex_x<<", " <<Jex_y<<", " <<Jex_z<<" )"<<std::endl;
	}	
	Jex_x*= 1e-3;  Jex_y*= 1e-3; Jex_z*= 1e-3; //pass to eV

	

	UnitCellK_TBOp //We usually need a set of operators 
	Ham(orbsPerCell),		//Hamiltonian
	VelX(orbsPerCell), 		//Velocity in X direction
	VelY(orbsPerCell), 		//Velocity in Y direction
	JXSz(orbsPerCell),		//Spin Current Operator in X direction
	JYSz(orbsPerCell),		//Spin Current Operator in Y direction
	Sx(orbsPerCell),			//Spin operator in x,y,z directions
	Sy(orbsPerCell),
	Sz(orbsPerCell),
	SxAB(orbsPerCell),		//Sublattice resolve Spin operator in x,y,z directions
	SyAB(orbsPerCell),
	SzAB(orbsPerCell);


	for(qt::integer io=0;io<num_orb ; io++)
	for(qt::integer is=0;is<max_spin; is++)
	{
		const qt::integer
		jo= 1-io,	//pseudospin flip index
		js= 1-is;	//spin flip index
			
		const qt::real
		oo= 1 - 2*io,//pseudospin sign
		ss= 1 - 2*is;//spin sign


		qt::integer
		row= io*max_spin + is,//variable for the rows
		col= row; 						//variable for the columns

		//define the indexes for 
		//the nearest and next-nearest neighbors indexes
		qt::integer 
			DI0nn[3]={0,-oo,0 },		
			DI1nn[3]={0,0 ,-oo},
			DI0nnn[6]={-1,+1, 0, 0 ,-1 ,+1 },
			DI1nnn[6]={ 0, 0,-1,+1 ,+1 ,-1 };
	
	//----------------------------ONSITE TERMS/----------------------------/
		{	//Include all the onsite terms
			col= io*max_spin + is;
			
			qt::complex 
			sz_val= ss ,
			h_val =
						StdPot*oo+		//scalar
						Jex_z*sz_val;	//Zeeman in Z

			Ham.AddEntry(row,col,h_val,0,0,0);
			Sz.AddEntry(row,col, ss , 0,0,0);
			SzAB.AddEntry(row,col, oo*ss, 0,0,0);
		}
	//----------------------------Lattice On-site with spin-flip/----------------------------/
		if(max_spin==2){
			col= io*max_spin + js;

			qt::complex
			sx_val= 1.0 ,
			sy_val=qt::complex(0.0,-ss),
			h_val = Jex_x*sx_val + Jex_y*sy_val;

			Ham.AddEntry(row,col,h_val,0,0,0);
			Sx.AddEntry(row,col,sx_val,0,0,0);
			Sy.AddEntry(row,col,sy_val,0,0,0);
			SxAB.AddEntry(row,col,oo*sx_val,0,0,0);
			SyAB.AddEntry(row,col,oo*sy_val,0,0,0);	

		}
	
	//----------------------------NEAREST NEIGHBORS no spinflip/----------------------------/
		for(qt::integer n=0;n<3;n++)
		{	
			const qt::real
			drsc[2]={ 
						DI0nn[n]*Lat[0][0] + DI1nn[n]*Lat[1][0] , 
						DI0nn[n]*Lat[0][1] + DI1nn[n]*Lat[1][1]
					},
			dr[2]={ 	
					drsc[0] + Opos[jo][0] - Opos[io][0]  ,
					drsc[1] + Opos[jo][1] - Opos[io][1] 
					};
		//----------------------------No spinflip/----------------------------/
			col= jo*max_spin + is;
			{	const qt::complex
				h_val  = hop,	//nearest neighbo hopping
				vx_val =-h_val*qt::complex(0.0,dr[0]),
				vy_val =-h_val*qt::complex(0.0,dr[1]) ;
				Ham.AddEntry( row,col, h_val ,DI0nn[n],DI1nn[n],0);
				VelX.AddEntry(row,col,vx_val ,DI0nn[n],DI1nn[n],0);
				VelY.AddEntry(row,col,vy_val ,DI0nn[n],DI1nn[n],0);
				if(max_spin==2)
				{
					JXSz.AddEntry(row,col,ss*vx_val,DI0nn[n],DI1nn[n],0);
					JYSz.AddEntry(row,col,ss*vy_val,DI0nn[n],DI1nn[n],0);				
				}
			}
	
			//----------------------------spinflip/----------------------------/
			col= jo*max_spin + js;
			if(max_spin==2)
			{	const qt::complex
				h_val = CrossProductDotZ(dr, is,js )*qt::complex(0.0,RSO*2./3.), //rashba 
				vx_val=-h_val*qt::complex(0.0,dr[0]),
				vy_val=-h_val*qt::complex(0.0,dr[1]);
				Ham.AddEntry( row,col, h_val,DI0nn[n],DI1nn[n],0);
				VelX.AddEntry(row,col,vx_val,DI0nn[n],DI1nn[n],0);
				VelY.AddEntry(row,col,vy_val,DI0nn[n],DI1nn[n],0);
			}
		}

	//----------------------------NEXT NEAREST NEIGHBORS/----------------------------/			
		const qt::real 
		t2= 1.0,
		VI[2]  ={ ISOA, ISOB },
		VPIA[2]={ PIAA, PIAB };	
		
		const qt::complex
		iso_val[6]={
				(t2+ qt::complex(0.0,ss*oo*VI[io]/3./sqrt(3.))),
				(t2- qt::complex(0.0,ss*oo*VI[io]/3./sqrt(3.))),
				(t2- qt::complex(0.0,ss*oo*VI[io]/3./sqrt(3.))),
				(t2+ qt::complex(0.0,ss*oo*VI[io]/3./sqrt(3.))),
				(t2- qt::complex(0.0,ss*oo*VI[io]/3./sqrt(3.))),
				(t2+ qt::complex(0.0,ss*oo*VI[io]/3./sqrt(3.)))
			   };
		for(qt::integer n=0;n<6;n++)
		{	
			const qt::real 
			dr[2]={ 
					 DI0nnn[n]*Lat[0][0] + DI1nnn[n]*Lat[1][0],
					 DI0nnn[n]*Lat[0][1] + DI1nnn[n]*Lat[1][1] 
					};

			// No spin flip
			col= io*max_spin + is;
			{	const qt::complex 
				h_val  = iso_val[n],
				vx_val =-h_val*qt::complex(0.0,dr[0]),
				vy_val =-h_val*qt::complex(0.0,dr[1]);
				Ham.AddEntry(	row,col,h_val ,DI0nnn[n],DI1nnn[n],0);
				VelX.AddEntry(	row,col,vx_val,DI0nnn[n],DI1nnn[n],0);
				VelY.AddEntry(	row,col,vy_val,DI0nnn[n],DI1nnn[n],0);
				if(max_spin==2)
				{
					JXSz.AddEntry(	row,col,ss*vx_val,DI0nnn[n],DI1nnn[n],0);
					JYSz.AddEntry(	row,col,ss*vy_val,DI0nnn[n],DI1nnn[n],0);
				}				
			}

			// Spin flip
			col= io*max_spin + js;
			if(max_spin==2)
			{	const qt::complex
				h_val = CrossProductDotZ(dr, is,js )*qt::complex(0.0,VPIA[io]*2./3.), //PIA
				vx_val=-h_val*dr[0] ,
				vy_val=-h_val*dr[1] ;
				VelX.AddEntry(row,col,vx_val,DI0nnn[n],DI1nnn[n],0);
				VelY.AddEntry(row,col,vy_val,DI0nnn[n],DI1nnn[n],0);			 
			}
		}
	}		

	Ham.writeOperator(label+".Ham");
	VelX.writeOperator(label+".VelX");
	VelY.writeOperator(label+".VelY");

	
	if(max_spin==2)
	{
		Sx.writeOperator(label+".Sx");
		Sy.writeOperator(label+".Sy");
		Sz.writeOperator(label+".Sz");
		SxAB.writeOperator(label+".SxAB");
		SyAB.writeOperator(label+".SyAB");
		SzAB.writeOperator(label+".SzAB");
		JXSz.writeOperator(label+".JxSz");
		JYSz.writeOperator(label+".JySz");
	}
	

return 0;}



