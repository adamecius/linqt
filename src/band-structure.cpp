// C libraries
// C++ libraries
// C libraries
// C++ libraries
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

// custom libraries
#include "types_definitions.hpp"
#include "parser.hpp"
#include "uck_tb_operator.hpp" // class SCTBOp
#include "dense_matrix.hpp"

//external libraries
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h"


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

	
	//Get Label
	std::string label, suffix;
	parser::GetValue(configFile,"Label" ,label );
	parser::GetValue(configFile,"suffix" ,suffix, parser::Optional); label+=suffix;

	//Get the number of path Kpath
	qt::integer kpathNum;
	parser::GetValue(configFile,"NumbeOfKPaths" ,kpathNum );

	std::vector< std::vector<qt::real> > kPath;
	parser::GetBlock(configFile, "KPATH", kPath , kpathNum ,4  );

	//Read the Operators labels
	std::string header = "OPERATORS";
	parser::FindHeader(configFile, header );

	std::string Ham;
	qt::integer max_spin;
	parser::GetValue(configFile,"Hamiltonian" , Ham );

	//Closing the input file
	configFile.close();

	UnitCellK_TBOp
	Hk( "operators/"+label+"."+Ham+".KOP" );

	
	std::cout<<"Computing the band structure of a system with lattice vectors: "<<std::endl;
	std::cout<<"LAT_VECTORS"<<std::endl;
	std::cout<<Hk.Lat(0,0)/Hk.LatConst()  <<" "<<Hk.Lat(0,1)/Hk.LatConst()  <<" "<<Hk.Lat(0,2)/Hk.LatConst()<<std::endl;
	std::cout<<Hk.Lat(1,0)/Hk.LatConst()  <<" "<<Hk.Lat(1,1)/Hk.LatConst()  <<" "<<Hk.Lat(1,2)/Hk.LatConst()<<std::endl;
	std::cout<<Hk.Lat(2,0)/Hk.LatConst()  <<" "<<Hk.Lat(2,1)/Hk.LatConst()  <<" "<<Hk.Lat(2,2)/Hk.LatConst()<<std::endl;
	std::cout<<"\nREC_VECTORS"<<std::endl;
	std::cout<<Hk.Rec(0,0)*Hk.LatConst()  <<" "<<Hk.Rec(0,1)*Hk.LatConst()  <<" "<<Hk.Rec(0,2)*Hk.LatConst()<<std::endl;
	std::cout<<Hk.Rec(1,0)*Hk.LatConst()  <<" "<<Hk.Rec(1,1)*Hk.LatConst()  <<" "<<Hk.Rec(1,2)*Hk.LatConst()<<std::endl;
	std::cout<<Hk.Rec(2,0)*Hk.LatConst()  <<" "<<Hk.Rec(2,1)*Hk.LatConst()  <<" "<<Hk.Rec(2,2)*Hk.LatConst()<<std::endl;
	
//Convert the path to cartesian units using the reciprocal lattice vectors
	for( qt::integer n=0; n< kpathNum; n++)
	{
		qt::real tmp[3]={ kPath[n][1], kPath[n][2], kPath[n][3] };
		for( int i=0; i< 3; i++)// the 3 is the number of spatial dimensions
		{
			kPath[n][1+i]=0;
			for( int j=0; j< 3; j++)// the 3 is the number of spatial dimensions
				kPath[n][1+i]+= tmp[j]*Hk.Rec(j,i);
		}
	}
	
	std::cout<<"\nTHROUGH "<< kpathNum<<" PATHS"<<std::endl;
	for(qt::integer n=0 ; n < kpathNum-1; n++)
		std::cout	<<"From:  "	<<kPath[n+0][1]<<" "<<kPath[n+0][2]	<<" "<<kPath[n+0][3] 
					<<" ----> "	<<kPath[n+1][1]<<" "<<kPath[n+1][2]	<<" "<<kPath[n+1][3]
					<<std::endl<<"Using: "<<kPath[n+1][0]<<" points."<<std::endl;

	//Structures needed for the exact diagonalization
	const qt::integer orbPerCell=Hk.Dim();
	qt::dense_matrix<qt::complex> 
		Eig(orbPerCell,orbPerCell),HkMat(orbPerCell,orbPerCell);	
	
	std::vector<qt::real> 
		Ek(orbPerCell);

	//compute the band structure
	std::ofstream output_file ( ( label+".band").c_str());
	qt::real norm =0 ;
	qt::real k0[3]={kPath[0][1] ,kPath[0][2],kPath[0][3]};
	for(qt::index n=0; n < kpathNum-1 ; n++ )
	{		
		qt::real K0[3] ={ kPath[n+0][1], kPath[n+0][2] , kPath[n+0][3] };			
		qt::real KF[3] ={ kPath[n+1][1], kPath[n+1][2] , kPath[n+1][3] };	
		
		const qt::integer np = kPath[n+1][0];
		for(qt::index p=0; p < np; p++ )
		{	
			qt::real t = (qt::real) p /(qt::real) np ;
			qt::real k[3] ={(1.0-t)*K0[0] + t*KF[0],
							(1.0-t)*K0[1] + t*KF[1],
							(1.0-t)*K0[2] + t*KF[2] 
						    };		
			Hk.Addk_dependence(k);
			Eig=Hk.Mat();

			Eig.Print();
			LAPACKE_zheev( LAPACK_ROW_MAJOR, 'N', 'U', orbPerCell,  Eig.beginPtr(), orbPerCell, &Ek[0] );	

			qt::real tmp_norm=0;
			for( qt::index i=0; i<3; i++)
			{
				tmp_norm +=	pow( k[i]-k0[i],2);
				k0[i]=k[i];
			}
			norm+=sqrt(tmp_norm);

			output_file<<norm*Hk.LatConst()<<" ";
			for(qt::index o =0; o< orbPerCell ; o++)
				output_file<<Ek[o]<<" ";
			output_file<<std::endl;
		}
	}
	output_file.close();

return 0;}




