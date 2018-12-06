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


qt::complex MatrixElem(	const qt::integer dim,
						const qt::dense_matrix<qt::complex> & A, 
						const qt::complex* PsiL,
						const qt::complex* PsiR)
{
	qt::complex output=0;; 
	
	for( qt::integer i=0;  i< dim ; i++)
	for( qt::integer j=0;  j< dim ; j++)
		output+= std::conj(PsiL[i])*A(i,j)*PsiR[i] ;
	return output;
}

qt::real MeanValue(const qt::integer dim,
					const qt::dense_matrix<qt::complex> & A,
					const qt::complex* Psi)
{
	qt::real output=0;; 
	
	for( qt::integer i=0;  i< dim ; i++)
	for( qt::integer j=0;  j< dim ; j++)
		output+= (std::conj(Psi[j])*A(j,i)*Psi[i]).real() ;
	return output;
}


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
	qt::integer TotPoint =0;
	qt::integer kpathNum;
	parser::GetValue(configFile,"NumbeOfKPaths" ,kpathNum );

	std::vector< std::vector<qt::real> > kPath;
	parser::GetBlock(configFile, "KPATH", kPath , kpathNum ,4  );

	for(int i=0; i< kpathNum; i++)
		TotPoint += kPath[i][0];
	
	std::cout<<"AQUI"<<TotPoint<<std::endl;

	std::string Ham;
	std::string ProjVec;
	std::vector<std::string> Proj;

	//Read the Operators labels
	std::string header = "OPERATORS";
	parser::FindHeader(configFile, header );
	parser::GetValue(configFile,"Hamiltonian" , Ham );
	parser::GetValue(configFile,"Projections", ProjVec);

	//tokenize string
	while( ProjVec.size() != 0 )
	{
		std::size_t
		pos = ProjVec.find(",");
		Proj.push_back( ProjVec.substr(0,pos) );
		if( pos==std::string::npos)
			ProjVec="";
		else 
			ProjVec = ProjVec.substr (pos+1,-1);
	} 
	qt::integer numproj = Proj.size();
	//Closing the input file
	configFile.close();

	UnitCellK_TBOp
	Hk( "operators/"+label+"."+Ham+".KOP" );
	std::vector< UnitCellK_TBOp > Ok( numproj );
	
	for(int i=0;i< numproj ; i++)
		Ok[i]= UnitCellK_TBOp( "operators/"+label+"."+Ham+".KOP" );

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

	std::vector< qt::real > k_norm(TotPoint+1,0);
	std::vector< qt::real > bands( TotPoint*(numproj+1)*orbPerCell);
	qt::real 
	k[3]= {kPath[0][1] ,kPath[0][2],kPath[0][3]};// initial point of kpath
		
	//In this loop we go through the selected K-Path


	//first one compute the initial point n==0
	qt::integer pidx = 0;
	{
		Hk.Addk_dependence(k);
		Eig=Hk.Mat();
		LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'U', orbPerCell,  Eig.beginPtr(), orbPerCell, &Ek[0] );	
		Eig.Transpose();
		//Loop for creating the band structure
		for(qt::index b =0; b< orbPerCell ; b++) //band index loop
		{
			//Add the energies to the band structure
			bands[( b*TotPoint + pidx)*(numproj+1)+0] = Ek[b];
			//Add the projections to the band structure
			for(qt::index p =0; p< numproj ; p++)
				bands[( b*TotPoint + pidx)*(numproj+1)+p+1] =  MeanValue(orbPerCell, Hk.Mat(),&Eig(b,0) );
		}
		pidx++;
	}
	//Then you compute the rest
	for(qt::index  n=1; n < kpathNum; n++ ) //Path loop
	for(qt::index lp=0; lp < kPath[n][0]; lp++ ) // local k-point for kpath
	{
		const qt::integer 
		np  = kPath[n][0]; 	//number of kpoint for a given path


		qt::real t = (qt::real)( lp+1 )/(qt::real)np;

		//We compute the norm
		for( int i=0;i<3 ;i++)
			k_norm[pidx]+= pow( ( kPath[n+1][i+1]-k[i] )*t, 2.0 );
		k_norm[pidx]= k_norm[pidx-1] + sqrt(k_norm[pidx]);

		//And  update both k and k0
		for( int i=0;i<3 && n < kpathNum;i++)
			k[i] = (1.0-t)*kPath[n][i+1] + t*kPath[n+1][i+1];


		Hk.Addk_dependence(k);
		Eig=Hk.Mat();
		LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'U', orbPerCell,  Eig.beginPtr(), orbPerCell, &Ek[0] );	
		Eig.Transpose();
		//Loop for creating the band structure
		for(qt::index b =0; b< orbPerCell ; b++) //band index loop
		{
			//Add the energies to the band structure
			bands[( b*TotPoint + pidx)*(numproj+1)+0] = Ek[b];
			//Add the projections to the band structure
			for(qt::index p =0; p< numproj ; p++)
				bands[( b*TotPoint + pidx)*(numproj+1)+p+1] =  MeanValue(orbPerCell, Hk.Mat(),&Eig(b,0) );
		}


		pidx++; 			//increace path point id to zero		
	}
	
	//Open the output file and print a header
	std::ofstream output_file ( ( label+".proj_band").c_str());
	output_file.precision(17);

	
	
	
	output_file<<"#  band_num    = "<<orbPerCell<<std::endl;
	output_file<<"#  NumKP    = "<<TotPoint<<std::endl;
	output_file<<"#  Projections = ";
	for(int i=0;i< numproj ; i++)
	output_file<<Proj[i]<<" ";
	output_file<<std::endl;

	//Now we print the bands
	for(qt::index b=0; b < orbPerCell ; b++ )
	{
		for(qt::index pidx=0; pidx < TotPoint 	; pidx++ )
		{
			output_file<<k_norm[pidx]<<" ";
			for(qt::index p=0; p < (1+numproj) 	; p++ )
				output_file<<bands[( b*TotPoint + pidx)*(1+numproj)+p]<<" ";
			output_file<<std::endl;
		}
		output_file<<std::endl;
	}		 
	output_file.close();

return 0;}




