
#include <mpi.h>
	// C libraries
// C++ libraries
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <iterator>


// custom libraries
#include "types_definitions.hpp"
#include "parser.hpp"
#include "kpm_dense_matrix.hpp"
#include "kspace-tboperator.h"


//external libraries
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h"


//The matrix convenction is
// AdwAdw AdwBdw AdwAup AdwBup
// BdwAdw BdwBdw BdwAup BdwBup
// AupAdw AupBdw AupAup AupBup
// BupAdw BupBdw BupAup BupBup
// the index convection is ( js*norb+ jo )		

int main(int argc, char *argv[])
{
//<<<<<<<<<<<<<<<< MPI INIT BLOCK >>>>>>>>>>>>>>>>>>>>>>>//
	MPI_Init(&argc,&argv);
	int world_rank, world_size, world_root=0;
	MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &world_size);		
	if( world_rank == world_root )
			std::cout<<"The program is running on "<<world_size<<" processes."<<std::endl; 
//<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
	if( argc != 2 )
	{
		std::cerr<<std::endl
				<<"Please submit a config input "<<std::endl;
		MPI_Finalize ( );
		return 0;
	}


	//If no exception, open the file 
	std::ifstream configFile;
	configFile.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
	try  {	configFile.open ( argv[1] );	}
	catch (std::ifstream::failure e)
	{
		configFile.close();
		std::cerr 	<<"Exception:"<<e.what()
					<<"\n while opening config file: "
					<< argv[1] <<"."<<std::endl;
		MPI_Finalize ( );
		return 0;
	}

	std::string label;
	GetValue(configFile,"Label" ,label );
	
	std::vector<int> kdim;
	GetBlock(configFile, "KPOINT", kdim );

	std::cout<<"Computing the density of statets for a grid "<<kdim[0]<<"x"<<kdim[1]<<"x"<<kdim[2]<<std::endl;

	//Closing the input file
	configFile.close();

	//INITALIZE THE OPERATORS
	std::string HAM_NAME="operators/"+label+".Ham.KOP";
	kspace::TBOperator HAM( HAM_NAME );
	HAM.readLatticeInfo( "operators/"+label+".LAT");	

	double b[3][3] ;
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		b[i][j]= HAM.Rec(i,j)/(double)kdim[i];

	int NE=4000;
	kpm::real Emin=-4.5,Emax=-Emin;	
	kpm::real En[NE],dos[NE];
	
	for(int n =0; n<NE; n++)
	{
		En[n] = Emin+ (double)( (double)n/(double)(NE-1) )*(Emax-Emin) ;
		dos[n]=0;
	}
	double eta=1./1000.;

	//Structures needed for the exact diagonalization
	int orbPerCell=HAM.Dimension();
	kpm::complex_matrix Eig(orbPerCell,orbPerCell),HamMat(orbPerCell,orbPerCell);	
	std::vector<kpm::real> Ek(orbPerCell);
	kpm::complex I = kpm::complex(0,1.0);


	size_t NumKP  = kdim[0]*kdim[1]*kdim[2] ;
	size_t myNumKP= ((long)NumKP + (long)world_size-1)/ (long)world_size ;
	size_t ko= world_rank*myNumKP ;
	size_t kf= ko +myNumKP ; if( kf>myNumKP	 ) kf= myNumKP ;

	for(size_t ik0=0; ik0 < kdim[0] ; ik0++ )
	for(size_t ik1=0; ik1 < kdim[1] ; ik1++ )	
	for(size_t ik2=0; ik2 < kdim[2] ; ik2++ )	
	{	
		const double sIdx = ( ik0*kdim[1] + ik1 )*kdim[2] + ik2 ; 
		if( sIdx >= ko && sIdx < kf ) 
		{
			const kpm::real
			k[3]    ={	ik0*b[0][0] + ik1*b[1][0] + ik1*b[2][0], 
						ik0*b[0][1] + ik1*b[1][1] + ik1*b[2][1],
						ik0*b[0][2] + ik1*b[1][2] + ik1*b[2][2]
					 } ;

			HAM.Addk_dependence(k);
			Eig=HAM.Mat();
			LAPACKE_zheev( LAPACK_ROW_MAJOR, 'N', 'U', orbPerCell,  Eig.beginPtr(), orbPerCell, &Ek[0] );	

			for(size_t n=0;n<NE;n++)
			for(size_t ei=0;ei<orbPerCell;ei++)
				dos[n] +=- ( 1.0/( En[n]- Ek[ei] +I* eta ) ).imag()/M_PI;
		}
	}


	kpm::real dosF[NE];
        MPI_Reduce( &dos[0], &dosF[0], NE, MPI_DOUBLE, MPI_SUM, world_root, MPI_COMM_WORLD);

	if( world_rank == world_root )
	{
		std::ofstream output_file("dos_test");
		for(size_t n=0;n<NE;n++)
			output_file<<En[n]<<" "<<dosF[n]/kdim[0]/kdim[1]<<std::endl;
		output_file.out;
	}

	MPI_Finalize ( );
return 0;}




