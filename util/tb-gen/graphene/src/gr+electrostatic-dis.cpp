
#include <cstdlib>
#include <sys/time.h>
#include <vector>
#include <algorithm>    // std::sort
#include <fstream>
#include <string>
#include <complex>
#include <iostream>
#include <iomanip>
#include "string_util.h"

#include "radial_disorder.hpp"
#include "lattice_geometry.h"
#include "mat_operators.hpp"
#include "tb_operator.hpp"
#include "MKS_physical_constants.hpp"

//SHARED FILES
#include "parser.hpp"
#include "validation.hpp"


static const Real tol_zero= std::numeric_limits<Real>::epsilon() ;
static const Real epsilon= std::numeric_limits<Real>::epsilon() ;

int IndexesToIndex( const int i0,const int  n0,
					const int i1,const int  n1,
					const int i2,const int  n2,
					const int i3,const int  n3 
					)
{
	return ( ( (i0+n0)%n0 *n1 + (i1+n1)%n1 )*n2 + (i2+n2)%n2 )*n3 + (i3+n3)%n3;  
}; // the index convection is ( jo*nspin + js )

void PrintProfile( 	const int  dim0, const int  dim1, const int  norb, const int spin,
					std::vector<double> Pot,
					std::string outputname)
{
	std::string SPINIDX[2]; SPINIDX[0]=".UP.DIS"; SPINIDX[1]=".DOWN.DIS";
	double MTONM=1E9;
	for(int is=0;is<spin ; is++)
	{
		std::ofstream outfile( (outputname+SPINIDX[is]).c_str() );
		for(int i0=0;i0<dim0   ; i0++)	
		for(int i1=0;i1<dim1   ; i1++)
		for(int io=0;io<norb ; io++)
		{	
			const int
			idx= IndexesToIndex(i0,dim0,i1,dim1,io,norb,is,spin);
			const Real 
			ri[3]={ i0*A[0][0]+ i1*A[1][0]+ io*Delta[0],
					i0*A[0][1]+ i1*A[1][1]+ io*Delta[1], 
					0
				  };
			outfile<<ri[0]*MTONM<<" "<<ri[1]*MTONM<<" "<<Pot[idx]<<std::endl;
		}		
		outfile.close();	
	}
};


int main(int argc, char *argv[])
{
	
	
	if( argc != 2)
	{
		std::cout<<"Please submit: label. The puddle, soc and exchange will be read out of label.inp "<<std::endl;
		return 0;
	}
	
        std::string header;
        std::string input_filename(argv[1]); 
        std::string label = input_filename.substr( 0,input_filename.rfind(".") ) ;
        std::cout<<label<<std::endl;
	//Check if the required files are present
	if( !validation::CheckFiles(1,&input_filename) )
	{
		std::cerr<<"ABORTING"<<std::endl;
		return -1;
	};
	std::ifstream configFile( input_filename.c_str() );


	
	int dim0=0, dim1=0, num_orb=0, max_spin=0;
	Real hop=0, StdPot;
	header="LATTICE_INFO";
	GetValue(configFile,"Dim0", dim0);
	GetValue(configFile,"Dim1", dim1);
	GetValue(configFile,"NumberOfOrbitals", num_orb);
	GetValue(configFile,"MaxSpin", max_spin);
	GetValue(configFile,"HoppingEnergy", hop);
	GetValue(configFile,"StaggeredPot", StdPot);
	
	
	std::cout	<<"Creating the hamiltonian for a lattice with:"<<std::endl
				<<" Dimensions:\t" << dim0<<"x"<<dim1<<std::endl
				<<" With NumberOfOrbitals:\t" << num_orb<<std::endl
				<<" and maximum spin:\t"<<max_spin<<std::endl
				<<" The leading hopping energy is :\t"<<hop<<" eV "<<std::endl
				<<" The staggered potential is :\t"<<StdPot<<" eV "<<std::endl;


	Real RSO, ISOA, ISOB, PIAA, PIAB;
	std::string soc_label;
	header="SOC_INFO";
	GetValue(configFile,"Spin-orbitName", soc_label);
        GetValue(configFile,"Rashba", RSO);
        GetValue(configFile,"IntrinsicA", ISOA);
        GetValue(configFile,"IntrinsicB", ISOB);
        GetValue(configFile,"PIA_A", PIAA);
        GetValue(configFile,"PIA_B", PIAB);

	Real Jex, Jtheta, Jphi, AF;
	header="MAGNET_INFO";
        GetValue(configFile,"Jex", Jex);
        GetValue(configFile,"J_Theta", Jtheta);
        GetValue(configFile,"J_Phi", Jphi);
        GetValue(configFile,"AF", AF);
		
	Real Jex_x= 0.001*Jex*sin(Jtheta*M_PI)*cos(Jphi*M_PI);
	Real Jex_y= 0.001*Jex*sin(Jtheta*M_PI)*sin(Jphi*M_PI);
	Real Jex_z= 0.001*Jex*cos(Jtheta*M_PI);

	const int dim= dim0*dim1*num_orb*max_spin;
	const int coord_max= 20;
	const Real TotArea=UnitCellArea*dim0*dim1; // en m
	const Real Ex = 0;

	//<<<<<<<<<<<<<<<<<<<<READ HAMILTONIAN PARAMETERS>>>>>>>>>>>>>>>>>>>>>>>>>//
	const Complex t1      = -hop;
	const Complex t2      =  0.0 ;
	const Complex VR      =  RSO/1000. ;
	const Complex VI[2]   =  { ISOA/1000.  , ISOB/1000.  };
	const Complex VPIA[2] =  { PIAA/1000.  , PIAB/1000.  };
	
	//<<<<<<<<<<<<<<<<<<<<ENDREAD HAMILTONIAN PARAMETERS>>>>>>>>>>>>>>>>>>>>>>//
	
	std::cout<<"Rashba    = "<< VR.real() <<std::endl;
	std::cout<<"Intrinsic = "<< ( VI[0].real() +  VI[1].real() )/2.0 <<std::endl;
	std::cout<<"VZ = "<< ( VI[0].real() -  VI[1].real() )/2.0 <<std::endl;
	std::cout<<"PIA = "<< VPIA[0].real() <<" "<< VPIA[1].real()<<std::endl;
	std::cout<<"Jex= ("<< Jex_x <<" , "<< Jex_y<<" , "<<Jex_z<<" ) "<<std::endl;
	std::cout<<"AF= "<<AF<<std::endl;


	TBOperator Ham(dim,coord_max);
	TBOperator VelX(dim,coord_max);
	TBOperator VelY(dim,coord_max);
	TBOperator VelXSz(dim,coord_max);
	TBOperator VelYSz(dim,coord_max);

	{
		TBOperator Sx(dim,1);
		TBOperator Sy(dim,1);
		TBOperator Sz(dim,1);
		TBOperator SxAB(dim,1);
		TBOperator SyAB(dim,1);
		TBOperator SzAB(dim,1);
		
		for(int i0=0;i0<dim0   ; i0++)
		for(int i1=0;i1<dim1   ; i1++)
		for(int io=0;io<num_orb ; io++)
		for(int is=0;is<max_spin; is++)
		{

			const int
			js= 1-is,
			row= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,is,max_spin);
			
			const Real
			oo= 2*io - 1.0  ;
			const Real
			ss= 1 - 2*is;

			int col;
			
			col= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,is,max_spin);
			Sz.AddEntry(row,col,ss/TotArea);
			SzAB.AddEntry(row,col,oo*ss/TotArea);
			
			col= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,js,max_spin);
			Sx.AddEntry(row,col,(Real)1.0/TotArea);
			Sy.AddEntry(row,col,-I*ss/TotArea);
			SxAB.AddEntry(row,col,(Real)oo/TotArea);
			SyAB.AddEntry(row,col,-I*oo*ss/TotArea);
		}
		Sx.WriteIntoFile(label+".Sx");
		Sy.WriteIntoFile(label+".Sy");
		Sz.WriteIntoFile(label+".Sz");
		SxAB.WriteIntoFile(label+".SxAB");
		SyAB.WriteIntoFile(label+".SyAB");
		SzAB.WriteIntoFile(label+".SzAB");
	}
	
	for(int i0=0;i0<dim0   ; i0++)	
	for(int i1=0;i1<dim1   ; i1++)
	for(int io=0;io<num_orb ; io++)
	for(int is=0;is<max_spin; is++)
	{	
		//Compute the sign associate with each index
		const Real
		ss= 1 - 2*is,	//is=(0)1 implies Sz=(1)-1
		oo=-1 + 2*io;
		//compute the final orbital/spin
		const int
		jo=1-io,
		js=1-is;
		//define the indexes for the nearest neighbors
		int DI0nn[3]={0,-oo,0 };
		int DI1nn[3]={0,0 ,-oo};
		//define the indexes of the next nearest neighbors
		int DI0nnn[6]={-1,+1, 0, 0 ,-1 ,+1 };
		int DI1nnn[6]={ 0, 0,-1,+1 ,+1 ,-1 };
		//define the inital index
		const int
		row= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,is,max_spin);
		const Real 
		ri[3]={ i0*A[0][0]+ i1*A[1][0]+ io*Delta[0],
				i0*A[0][1]+ i1*A[1][1]+ io*Delta[1], 
				0
			  };
	//----------------------------ONSITE ENERGY/----------------------------/
		{
			const int
			col= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,is, max_spin);
			Complex val =
					StdPot*oo	+ //scalar
					Jex_z*ss*(1.0-2.0*io*AF); //Zeeman in Z
			Ham.AddEntry(row,col,val);
		}
	//----------------------------Spin-Flip ONSITE ENERGY/----------------------------/
		{
			const int
			col= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,js, max_spin);
			Complex val = 	Jex_x -
					Jex_y*ss*I;
			Ham.AddEntry(row,col,val);
		}
	//----------------------------NEAREST NEIGHBORS/----------------------------/
		for(int n=0;n<3;n++)
		{	
			const Real 
			rj[2]={ (i0+DI0nn[n])*A[0][0] + (i1+DI1nn[n])*A[1][0]+jo*Delta[0],
					(i0+DI0nn[n])*A[0][1] + (i1+DI1nn[n])*A[1][1]+jo*Delta[1]
					};
			const int
			col= IndexesToIndex(i0+DI0nn[n],dim0,i1+DI1nn[n],dim1,jo,num_orb,is,max_spin);
			Complex val = t1;
			Ham.AddEntry(row,col,val);
			VelX.AddEntry(row,col,-I*(rj[0]-ri[0])*val);
			VelY.AddEntry(row,col,-I*(rj[1]-ri[1])*val);
			VelXSz.AddEntry(row,col,-ss*I*(rj[0]-ri[0])*val);
			VelYSz.AddEntry(row,col,-ss*I*(rj[1]-ri[1])*val);				
		}
	//----------------------------NEXT NEAREST NEIGHBORS/----------------------------/			
		const Complex
		iso_val[6]={
				(t2+ I*ss*oo*VI[io]/3./sqrt(3.)),
				(t2- I*ss*oo*VI[io]/3./sqrt(3.)),
				(t2- I*ss*oo*VI[io]/3./sqrt(3.)),
				(t2+ I*ss*oo*VI[io]/3./sqrt(3.)),
				(t2- I*ss*oo*VI[io]/3./sqrt(3.)),
				(t2+ I*ss*oo*VI[io]/3./sqrt(3.))
			   };
		for(int n=0;n<6;n++)
		{	
			const Real 
			rj[2]={ 
					( i0+DI0nnn[n] )*A[0][0] + (i1+DI1nnn[n])*A[1][0] + io*Delta[0],
					( i0+DI0nnn[n] )*A[0][1] + (i1+DI1nnn[n])*A[1][1] + io*Delta[1] 
				};

			const int
			col= IndexesToIndex(i0+DI0nnn[n],dim0,i1+DI1nnn[n],dim1,io,num_orb,is,max_spin);
			Complex val = iso_val[n];
			Ham.AddEntry(row,col,val);
			VelX.AddEntry(row,col,-I*(rj[0]-ri[0])*val);
			VelY.AddEntry(row,col,-I*(rj[1]-ri[1])*val);
			VelXSz.AddEntry(row,col,-ss*I*(rj[0]-ri[0])*val);
			VelYSz.AddEntry(row,col,-ss*I*(rj[1]-ri[1])*val);
		}		
	//---------------------------- NEAREST NEIGHBORS SPIN FLIP/----------------------------/
		for(int n=0;n<3;++n)
		{
			
			//compute the r_j position
			const Real
			rj[2]=	{
					(i0+DI0nn[n])*A[0][0] + (i1+DI1nn[n])*A[1][0] + jo*Delta[0],
					(i0+DI0nn[n])*A[0][1] + (i1+DI1nn[n])*A[1][1] + jo*Delta[1]
				};
			//compute the differen between ri and rj
			const Real
			del_r[3]={ ri[0]-rj[0] ,ri[1]-rj[1] ,0.  };
			//compute the rashba value
			const Complex
			PhiRSO= 2.*I*VR*CrossProductDotZ(del_r, is,js )/3.;

			const int
			col= IndexesToIndex(i0+DI0nn[n],dim0,i1+DI1nn[n],dim1,jo,num_orb,js,max_spin);
			Complex val = PhiRSO;
			Ham.AddEntry(row,col,val);
			VelX.AddEntry(row,col,-I*(rj[0]-ri[0])*val);
			VelY.AddEntry(row,col,-I*(rj[1]-ri[1])*val);
		}
	//---------------------------- NEXT  NEAREST NEIGHBORS SPIN FLIP/-----------------------
			//------ PIA INTERACTION------
			for(int n=0;n<6;++n)
			{
				const Real 
				rj[2]={ ( i0+DI0nnn[n] )*A[0][0] + ( i1+DI1nnn[n] )*A[1][0]+ io*Delta[0],
						( i0+DI0nnn[n] )*A[0][1] + ( i1+DI1nnn[n] )*A[1][1]+ io*Delta[1] 
						};						
				//compute the differen between ri and rj					
				const Real 
				del_r[3]={ ri[0]-rj[0] ,ri[1]-rj[1] ,0.  };
				//Compute Pia phase
				const Complex 
				PhiPIA= 2.*I*VPIA[io]*CrossProductDotZ(del_r, is,js )/3.;

				const int
				col= IndexesToIndex(i0+DI0nnn[n],dim0,i1+DI1nnn[n],dim1,io,num_orb,js,max_spin);
				Complex val = PhiPIA;
				Ham.AddEntry(row,col,val);
				VelX.AddEntry(row,col,-I*(rj[0]-ri[0])*val);
				VelY.AddEntry(row,col,-I*(rj[1]-ri[1])*val);



			}		

	}

	//Add diagonal disorder
	std::vector<Real> diagElements(dim,0.);
	bool dis_found;

	std::cout<<std::endl<<"Adding disorder to the system:"<<std::endl;
		std::cout<<"-----------------------------------------------"<<std::endl<<std::endl;
	
	header = "ANDERSON_DISORDER";
	dis_found =(bool)( FindHeader(configFile,header) != string::npos );	
	if( dis_found )
	{
		Real dis_par[2]; int Seed;
		GetArray(configFile,"DisorderParameters", &dis_par[0],1,header,parser::NotOptional); 
		std::cout	<<"Disorder: "<<header<<" was found."<< std::endl
					<<"\tStrength: "<<dis_par[0]<<" eV"<< std::endl ;
		std::cout<<"-----------------------------------------------"<<std::endl<<std::endl;
	}

	header = "STAGGERED_DISKS";
	dis_found =(bool)( FindHeader(configFile,header) != string::npos );	
	if( dis_found )
	{
		Real dis_par[3]; int Seed=0;
		GetArray(configFile,"DisorderParameters", &dis_par[0],3,header,parser::NotOptional); 
		GetValue(configFile,"Seed", Seed,parser::CurrLine,parser::Optional); 
		std::cout	<<"Disorder: "<<header<<" was found."<< std::endl
					<<"\tRadius: "<<dis_par[0]<<" nm"<< std::endl 
					<<"\tConcentration: "<<dis_par[1]<<" \%"<< std::endl 
					<<"\tStrength: "<<dis_par[2]<<" eV"<< std::endl 
					<<"Using Seed: "<<Seed<<std::endl;
		std::cout<<"-----------------------------------------------"<<std::endl<<std::endl;

		dis_par[0]*=1.0E-9;
		dis_par[1]/=100;
		dis_par[2];

		std::vector<Real> tmp(dim,0.);
		StaggeredDiskDisorder stgdis(dis_par[2],dis_par[1],dis_par[0]);		
		radial_disorder(	dim0, dim1, num_orb, max_spin, 
							stgdis,tmp,Seed
						);
		PrintProfile(dim0,dim1,num_orb,max_spin,tmp, label+"-"+header );
		for(int i=0;i<dim;i++)
			diagElements[i]+= tmp[i];
	}
	
	header = "CHARGE_PUDDLES";
	dis_found =(bool)( FindHeader(configFile,header) != string::npos );	
	if( dis_found )
	{
		Real dis_par[3]; int Seed=0;
		GetArray(configFile,"DisorderParameters", &dis_par[0],3,header,parser::NotOptional); 
		GetValue(configFile,"Seed", Seed,parser::CurrLine,parser::Optional); 
		std::cout	<<"Disorder: "<<header<<" was found."<< std::endl
					<<"\tRadius: "<<dis_par[0]<<" nm"<< std::endl 
					<<"\tConcentration: "<<dis_par[1]<<" \%"<< std::endl 
					<<"\tStrength: "<<dis_par[2]<<" eV"<< std::endl 
					<<"Using Seed: "<<Seed<<std::endl;
		std::cout<<"-----------------------------------------------"<<std::endl<<std::endl;

		dis_par[0]*=1.0E-9;
		dis_par[1]/=100;
		dis_par[2];

		std::vector<Real> tmp(dim,0.);
		PuddleDisorder pudd(dis_par[2],dis_par[1],dis_par[0]);		
		radial_disorder(	dim0, dim1, num_orb, max_spin, 
							pudd,tmp,Seed
						);
		PrintProfile(dim0,dim1,num_orb,max_spin,tmp, label+"-"+header );
		for(int i=0;i<dim;i++)
			diagElements[i]+= tmp[i];
	}

	for(int i0=0;i0<dim0   ; i0++ )
	for(int i1=0;i1<dim1   ; i1++ )
	for(int io=0;io<num_orb ; io++ )
	for(int is=0;is<max_spin ; is++ )
	{	
		const int
		k0= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,is, max_spin);
		Complex val = diagElements[k0];
		Ham.AddEntry(k0,k0,val);
	}

	std::cout<<"Checking the hermiticity of the matrix"<<std::endl;
	bool is_hermitian=true;
	for(int i0=0;i0<dim0   ; i0++ )
	for(int i1=0;i1<dim1   ; i1++ )
	for(int io=0;io<num_orb ; io++ )
	for(int is=0;is<max_spin ; is++ )
	{
		const int 
		row= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,is,max_spin),
		coordnum0=Ham.Op_triplet[row].size();

		for(int n=0; n<coordnum0 ; n++)
		{
			const int 
			col= Ham.Op_triplet[row][n].col;
			const Complex 
			val= Ham.Op_triplet[row][n].val;

			const int 
			coordnum1=Ham.Op_triplet[col].size();
			for(int m=0; m<coordnum1 ; m++)
				if ( Ham.Op_triplet[col][m].col == row )
					if ( std::norm( std::conj(Ham.Op_triplet[col][m].val) -val) > tol_zero )
						std::cout<<row<<" "<<col<<" "<<val<<" | "<<col<<" "<<Ham.Op_triplet[col][m].col<<" "<<  Ham.Op_triplet[col][m].val<<std::endl;
 

		}

	}
	Ham.WriteIntoFile(label+".Ham");


{
	const kpm::real E0=1.0;	// V/m
	const kpm::real Q0=1.0;	// in units of the electrons charge
	const kpm::real hb=0.6582119; 	// eV.fs
	const kpm::real ConvFact= E0*Q0/hb;
	for(int i0=0;i0<dim0   ; i0++ ) for(int i1=0;i1<dim1   ; i1++ )
	for(int io=0;io<num_orb ; io++ ) for(int is=0;is<max_spin ; is++ )
	{
		const int 
		row		= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,is,max_spin),
		coordnum=VelX.Op_triplet[row].size();
		for(int n=0; n<coordnum ; n++)
		VelX.Op_triplet[row][n].val = VelX.Op_triplet[row][n].val *ConvFact;

	}
	VelX.WriteIntoFile(label+".Pow");
	for(int i0=0;i0<dim0   ; i0++ )
	for(int i1=0;i1<dim1   ; i1++ )
	for(int io=0;io<num_orb ; io++ )
	for(int is=0;is<max_spin ; is++ )
	{
		const int 
		row		= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,is,max_spin),
		coordnum=VelX.Op_triplet[row].size();
		for(int n=0; n<coordnum ; n++)
			VelX.Op_triplet[row][n].val = VelX.Op_triplet[row][n].val /ConvFact;
	}
}

	//The velocities are in units of  [H,X] = eV m , therefore we first divide by hbar
	//The velocities are in units of  eV m , therefore we first divide by hbar
	const kpm::real Q0=1.602177E-19;	// Coulomb
	const kpm::real hb=6.582119E-16; 	// eV.s
	const kpm::real ConvFact= Q0/(TotArea*hb);
	std::cout<<"THE AREA IS: "<<TotArea<<" m^"<<std::endl;
	std::cout<<"Conversion factor"<<ConvFact<<std::endl;
	for(int i0=0;i0<dim0   ; i0++ )
	for(int i1=0;i1<dim1   ; i1++ )
	for(int io=0;io<num_orb ; io++ )
	for(int is=0;is<max_spin ; is++ )
	{
		const int
		row= IndexesToIndex(i0,dim0,i1,dim1,io,num_orb,is,max_spin);

		int coordnum=VelX.Op_triplet[row].size();
		for(int n=0; n<coordnum ; n++)
			VelX.Op_triplet[row][n].val*=ConvFact;

		coordnum=VelY.Op_triplet[row].size();
		for(int n=0; n<coordnum ; n++)
			VelY.Op_triplet[row][n].val*=ConvFact;

		coordnum=VelXSz.Op_triplet[row].size();
		for(int n=0; n<coordnum ; n++)
			VelXSz.Op_triplet[row][n].val*=ConvFact;

		coordnum=VelYSz.Op_triplet[row].size();
		for(int n=0; n<coordnum ; n++)
			VelYSz.Op_triplet[row][n].val*=ConvFact;
	}
	VelX.WriteIntoFile(label+".Jx");
	VelY.WriteIntoFile(label+".Jy");
	VelXSz.WriteIntoFile(label+".JxSz");
	VelYSz.WriteIntoFile(label+".JySz");

return 0;}



