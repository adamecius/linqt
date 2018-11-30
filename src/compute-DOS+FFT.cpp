
// LOAD MPI LIBRARY IF AVAILABLE
#include <mpi.h> //This should be the first header ALWAYS
// C libraries
// C++ libraries
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <iterator>

// External libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h> 
#include <sys/types.h>
#include <unistd.h>

// custom libraries
#include "types_definitions.hpp"
#include "kpm_config.hpp"	// class kpm::config
#include "uck_tb_operator.hpp" // class SCTBOp
#include "lattice_fftw3_mpi.hpp" // class LatticeFFTW
#include "kpm_block_matrix.hpp"
#include "kpm_linalg.hpp"
#include "kpm.hpp"

int main(int argc, char *argv[])
{

//<<<<<<<<<<<<<<<< MPI INIT BLOCK >>>>>>>>>>>>>>>>>>>>>>>//
	bool sucess ;
	MPI_Init(&argc,&argv);
	int world_rank, world_size, world_root=0;
	MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &world_size);		
	if( world_rank == world_root )
			std::cout<<"The program is running on "<<world_size<<" processes."<<std::endl; 
//<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//

	//Set the system label and try to open the configuration file
	// If fail to open it, abort the program
	std::string configFilename; 
	if ( argc!= 1) configFilename= argv[1];
	std::ifstream configFile( configFilename.c_str() );
	try{ 
		configFile.exceptions(configFile.failbit);
	}catch (const std::ios_base::failure& e){
		if ( world_rank == world_root)
			std::cerr 	<<"Exception:"<<e.what()
						<<"\n while opening config file: "
						<< configFilename <<std::endl;
		MPI_Finalize();
		return 0;
	}
//----------------------START CONFIG READING---------------------//
	bool valid_opt;
	kpm::config opt; 

	valid_opt = kpm::GetOptionsFromCFG(configFile,opt)*kpm::ValidateOptions(opt);
	if( world_rank==world_root)
	{
		if( !valid_opt)
		std::cerr	<<"Problems reading the data from: "
						<<configFilename<<std::endl
						<<"Maybe one the data field has incorrect values"
						<<std::endl;
		PrintOptions(opt);
	}
	if( !valid_opt)
	{
		MPI_Finalize();
		return 0;
	}

	//Read K-POINT INFORMATION
	std::vector<int> kdim;
	GetBlock(configFile, "KPOINT", kdim );

	//Read SUPERCELL INFORMATION
	std::vector<int> scdim;
	GetBlock(configFile, "SUPERCELL", scdim );
	double fft_norm = 1.0/scdim[0]/scdim[1]/scdim[2];
	
	//Read the Operators labels
	std::string header = "OPERATORS";
	std::string Ham;
	FindHeader(configFile, header );
	GetValue(configFile,"Hamiltonian" , Ham );

	//Read the Operators labels
	header = "RANDOM_VECTOR";
	FindHeader(configFile, header );
	int RandNum, RSeed;
	GetValue(configFile,"Sampling" , RandNum );
	GetValue(configFile,"RandomSeed" , RSeed );

	configFile.close();
//----------------------FINISH CONFIG READING---------------------//
	
	//For convention we define the following variables
	kpm::real EnergyScalFact = opt.energyScale /2.0*opt.cutOff ;
	const int MomSize=opt.maxMom;

	//Read UnitCell k-dependent Tight-binding Operators. 
	// and their lattice information
	UnitCellK_TBOp
	Hk( "operators/"+opt.label+"."+Ham+".KOP" );
	Hk.readLatticeInfo( "operators/"+opt.label+".LAT");	
	const int orbPerCell = Hk.Rank();
	const long int numCell = scdim[0]*scdim[1]*scdim[2];

	// Here we construct the equivalent reciprocal space vector
	// for the particular unit cell we read on K_UnitCell
	// but assuming we expanding into a super cell
	double b[3][3] ;
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		b[i][j]= Hk.Rec(i,j)/(double)scdim[i];

	// Finally we print the options to check everything is ok
	if( world_rank == world_root)
	{
		std::cout<<std::endl;
		std::cout<<"OPERATORS"<<std::endl;
		std::cout<<"The total number of orbitals per cell is: "<<orbPerCell<<std::endl;
		std::cout<<"The Hamiltonian label is "<<Ham<<std::endl;
		std::cout<<std::endl<<std::endl;
		std::cout<<"The algorithm is going to compute: "
				 <<MomSize<<" Chebyshev moments for the Density of states"<<std::endl;
		std::cout<<"RANDOM_VECTOR"<<std::endl;
		std::cout<<"The sampling number is "<<RandNum<<std::endl;
		std::cout<<"The random seed is "<<RSeed<<std::endl;
		std::cout<<"GEOMETRICAL_INFO"<<std::endl;
		std::cout<<"NUMBER OF ORBITALS: "<<orbPerCell<<std::endl;
		std::cout<<"NUMBER OF KPOINTS: "<<kdim[0]<<" x "<<kdim[1]<<" x "<<kdim[2]<<std::endl;
		std::cout<<"SUPERCELL DIMENSIONS="<<scdim[0]<<" x "<<scdim[1]<<" x "<<scdim[2]<<std::endl;
		std::cout<<"LAT_VECTORS [ LatConst ] "<<std::endl;
		std::cout<<Hk.Lat(0,0)/Hk.LatConst()<<" "<<Hk.Lat(0,1)/Hk.LatConst()<<" "<<Hk.Lat(0,2)/Hk.LatConst()<<std::endl;
		std::cout<<Hk.Lat(1,0)/Hk.LatConst()<<" "<<Hk.Lat(1,1)/Hk.LatConst()<<" "<<Hk.Lat(1,2)/Hk.LatConst()<<std::endl;
		std::cout<<Hk.Lat(2,0)/Hk.LatConst()<<" "<<Hk.Lat(2,1)/Hk.LatConst()<<" "<<Hk.Lat(2,2)/Hk.LatConst()<<std::endl;
		std::cout<<"REC_LAT_VECTORS (2pi/LatConst )"<<std::endl;
		std::cout<<Hk.Rec(0,0)*0.5*Hk.LatConst()/M_PI<<" "<<Hk.Rec(0,1)*0.5*Hk.LatConst()/M_PI<<" "<<Hk.Rec(0,2)*0.5*Hk.LatConst()/M_PI<<std::endl;
		std::cout<<Hk.Rec(1,0)*0.5*Hk.LatConst()/M_PI<<" "<<Hk.Rec(1,1)*0.5*Hk.LatConst()/M_PI<<" "<<Hk.Rec(1,2)*0.5*Hk.LatConst()/M_PI<<std::endl;
		std::cout<<Hk.Rec(2,0)*0.5*Hk.LatConst()/M_PI<<" "<<Hk.Rec(2,1)*0.5*Hk.LatConst()/M_PI<<" "<<Hk.Rec(2,2)*0.5*Hk.LatConst()/M_PI<<std::endl;
	}
	// Rescal the operators with their appropiate geometric factors
	// The hamiltonian should be rescale for KPM
	Hk.Rescale(1.0/EnergyScalFact);

	// The class LatticeFFTW will be used 
	// to control the parallelization over k
	LatticeFFTW PotFFT;
	sucess = PotFFT.Init(MPI_COMM_WORLD,scdim, orbPerCell,opt.label);
	// Read the opetions from label
	if(	!sucess 	)
	{ MPI_Finalize(); return 0; }

	const int 
	loc_dim = PotFFT.LocalDim(), 	//Local dimension dimension of the subspace per process 
	loc_scdim0 = PotFFT.LocalDim0(),//Local supercell dimension from a0 lat-vect
	i0_shift = PotFFT.LocalI0();	//coordinate shift due to the parallelization

	//Before proceding, we should first used the seed for initialize
	//the variables for gsl random vector
	int par_seed = RSeed* getpid()*( world_rank + 2321);    //create a seed that is different for process.
	gsl_rng_env_setup();
	gsl_rng *rng;  // random number generator
	rng = gsl_rng_alloc (gsl_rng_default);     // uses the default
//	gsl_rng_set (rng, par_seed);                  // set seed
//	gsl_rng_set (rng, RSeed); //FOR TESTING
	//Allocate the memory for all the necessary arrays (in parallel)
	std::vector<kpm::complex> data( 5*loc_dim );
	kpm::complex
		*Psi = &data[0*loc_dim],
		*jm0 = &data[1*loc_dim],
		*jm1 = &data[2*loc_dim],
		*jmt = &data[3*loc_dim],
		*jmU = &data[4*loc_dim];

	//Allocate the moment vector which is of the same size for all processes
	std::vector<kpm::complex> mu( MomSize,0.0 );

	//Then we read the Potential Operator
	//We will assume that this operator is written also as 
	//a block diagonal operator but in real space,
	kpm::BlockMatrix<kpm::complex> Uloc( loc_dim , orbPerCell );

	//Add Electrostatic Disorder
//	kpm::BlockMatrix<kpm::complex> Uloc( loc_dim , orbPerCell );

	//Add Structural Disorder
//	kpm::BlockMatrix<kpm::complex> Uloc( loc_dim , orbPerCell );

	//Add Adatoms
//	kpm::BlockMatrix<kpm::complex> Uloc( loc_dim , orbPerCell );

	for(size_t loc_i0=0; loc_i0 < loc_scdim0 ; loc_i0++ )
	for(size_t i1=0; i1 < scdim[1] 	 ; i1++ )	
	for(size_t i2=0; i2 < scdim[2] 	 ; i2++ )	
	{
		const kpm::real W = 2.8;
		const size_t bidx = ( loc_i0*scdim[1] + i1 )*scdim[2] + i2;
		for( int o=0; o<orbPerCell; o++ )
			Uloc(bidx,o,o) = 0.5*W*( 2.0*gsl_rng_uniform(rng) -1.0 )/EnergyScalFact;   
	}
//	Uloc.Print();
	
	//In addition, we create three BlockDiag Operators 
	//for the hamiltonian, the force, and the operators
	kpm::BlockMatrix<kpm::complex>
	Hloc( loc_dim , orbPerCell );

	//There are two possible loops, one is the K-point look
	//asociated with the periodic system, and the other is the
	// Nsc-dimensional loop in the momentum space associated with the 
	// the super cell of dimension Nsc. We begin by considering the 
	// Periodic loop, for which we will use the the indexes
	// ik0, ik1, and ik2;
	for(size_t ik0=0; ik0 < kdim[0] ; ik0++ )
	for(size_t ik1=0; ik1 < kdim[1] ; ik1++ )	
	for(size_t ik2=0; ik2 < kdim[2] ; ik2++ )	
	{

		//For each  (ik0,ik1,ik2) Define a k-shift dk in the BZ
		//defined by the periodicity of the unit cell
		// given by:
 		const kpm::real
		db[3] = { ik0/(double)kdim[0], ik1/(double)kdim[1]+ ik2/(double)kdim[2] } ;

		if( world_rank == world_root)
			std::cout	<<"Computing the moments for the k-shift: "<<std::endl
						<<db[0]<<" "<<db[1]<<" "<<db[2]<<std::endl;


		// For the super-cell we will use the indexes
		// i0,i1, and i2 and due to the dimensionality of the supercell
		// each Delta KBZ, will be sub-divided according to the number
		// of sites in the super cell.
		// This will be taking into account 
		// when creatking the local H,O and F matrices
		//which are assume to be block diagonal and therefore can have
		//and special format
		if( world_rank == world_root)
			std::cout	<<"Creating a matrix with "<<loc_dim/orbPerCell<<" blocks "
						<<"with dimension: "<< orbPerCell<<"x"<<orbPerCell<<std::endl
						<<"at each of the "<<world_size<<" nodes."<<std::endl;
		for(size_t loc_i0=0; loc_i0 < loc_scdim0 ; loc_i0++ )
		for(size_t i1=0; i1 < scdim[1] 	 ; i1++ )	
		for(size_t i2=0; i2 < scdim[2] 	 ; i2++ )	
		{

			const size_t bidx = ( loc_i0*scdim[1] + i1 )*scdim[2] + i2;
			const kpm::real
			k[3] = {
					(loc_i0+i0_shift + db[0])*b[0][0] + (i1 + db[1])*b[1][0]+ (i2+ db[2])*b[2][0],
					(loc_i0+i0_shift + db[0])*b[0][1] + (i1 + db[1])*b[1][1]+ (i2+ db[2])*b[2][1],
					(loc_i0+i0_shift + db[0])*b[0][2] + (i1 + db[1])*b[1][2]+ (i2+ db[2])*b[2][2]
					};

			Hk.Addk_dependence(k );
			Hloc.AddBlock( bidx , Hk.Mat() ) ;
 
		}
	
		//Before Starting KPM we define the vector for the trace
		//We make no assumption on anything for the moment
		for( int r=0; r <  RandNum ; r++ )
		{
			// The random vector will varie from kBZ point
			// to kBZ point to simulate a larger vector
			for( int n=0 ; n< loc_dim ; n++)
			{
				double theta_r=  ( 2.0*gsl_rng_uniform(rng) -1.0 )*M_PI;
				Psi[ n ] = std::complex<double>( cos( theta_r ), sin( theta_r) );
			}

			//------------------------WE NOW START KPM--------------------/
			kpm::copy(loc_dim,jm0,Psi);
			for(int  m=0;m< MomSize ; m++ )
			{
				//Luego iteramos normalmente
				kpm::complex mutmp;
				// T_m(H) | POW Psi>
				cheb_evolve( m , Hloc,jm0, jm1 );

				//We compute the potential effect and reserve the 
				//vector
				kpm::real b=2; if ( m == 0 ) b=1;
				kpm::copy(loc_dim,jmt,jm0);
				
				PotFFT.DirectFourierTransform(jmt);
				Uloc.Multiply( b , jmt , 0 , jmU ); 
				PotFFT.InverseFourierTransform(jmU);
//				//Here we add the potential effect
				kpm::axpy(loc_dim,fft_norm, jmU,jm1);

				//Add the multiplication  with the potential
				kpm::dot(loc_dim, Psi , jm0 ,mutmp); //<Psi T_m(H) Op | T_m(H) Pow Psi >
				mu[ m ] += mutmp ;
			}
		}
	}

		
	//After we finished with all the moments, each cell has parts
	//of the final answer. We will collect all of these in cell 0
	std::vector<kpm::complex> mu_Final(MomSize,0.0);	
	MPI_Reduce( 
				&mu[0]		,  &mu_Final[0], 
				MomSize, MPI_DOUBLE_COMPLEX, 
				MPI_SUM, world_root, MPI_COMM_WORLD
				); 

	// Finally save all the moments
	if( world_rank == world_root )
	{
		std::ofstream mom_file(( opt.label+opt.suffix+".DOSmom").c_str());

		mom_file	<<opt.label+opt.suffix<<" "<<MomSize<<" " //Label  NumChebMoms
					<<opt.energyShift- opt.energyScale*0.5/opt.cutOff<<" " 	// Lower end of the band
					<<opt.energyShift+ opt.energyScale*0.5/opt.cutOff<<" " 	//Higher end of the band
					<<opt.cutOff<<" "<<1<<" "<<1<<" "<<1<<" "<<std::endl; //Options to be depreciated

		const double hbar= 0.6582119E-15;
		for( int m0=0;m0< MomSize ; m0++ )
		{
			kpm::complex muelem= mu_Final[ m0]/scdim[0]/scdim[1]/scdim[2];
			mom_file<<m0<<" "<<muelem.real()<<" "<<muelem.imag()<<std::endl;
		}
		mom_file.close();
	}


	MPI_Finalize();
return 0;}
