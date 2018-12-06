
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
#include "kpm.hpp"	// class kpm::config
#include "parser.hpp"
#include "uck_tb_operator.hpp" // class SCTBOp
#include "lattice_fftw3_mpi.hpp" // class LatticeFFTW


#include "sparse_block_matrix.hpp"
#include "vector_linalg.hpp"




int main(int argc, char *argv[])
{
//<<<<<<<<<<<<<<<< MPI INIT BLOCK >>>>>>>>>>>>>>>>>>>>>>>//
	
	int world_rank, world_size, world_root=0;
	MPI_Init(&argc,&argv);
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
	kpm::GetOptionsFromCFG(configFile,opt);
	
	if( world_rank == 0 )
		kpm::PrintOptions(opt);


	//Read K-POINT INFORMATION
	std::vector<qt::dimension> kdim;
	parser::GetVector(configFile, "KPOINT", kdim ,3	);


	//Read SUPERCELL INFORMATION
	std::vector<qt::dimension> scdim;
	parser::GetVector(configFile, "SUPERCELL", scdim ,3 );
	qt::real fft_norm = 1.0/scdim[0]/scdim[1]/scdim[2];


	//Read the Operators labels
	std::string header = "OPERATORS";
	std::string Ham, For, Obs;
	parser::FindHeader(configFile, header );
	parser::GetValue(configFile,"Hamiltonian" , Ham );
	parser::GetValue(configFile,"Force", For );
	parser::GetValue(configFile,"Observable", Obs );


	//Read the Operators labels
	header = "RANDOM_VECTOR";
	qt::integer RandNum, RSeed;
	parser::FindHeader(configFile, header );
	parser::GetValue(configFile,"Sampling" , RandNum );
	parser::GetValue(configFile,"RandomSeed" , RSeed );

	configFile.close();



//----------------------FINISH CONFIG READING---------------------//
	
	//For convention we define the following variables
	const qt::integer 
		MomSize = (qt::integer)opt.maxMom ;

	const qt::real 
		EnergyScalFact = opt.energyScale /2.0*opt.cutOff ;


	//Read UnitCell k-dependent Tight-binding Operators. 
	// and their lattice information

	UnitCellK_TBOp
	Hk( "operators/"+opt.label+"."+Ham+".KOP" ),
	Fk( "operators/"+opt.label+"."+For+".KOP" ),
	Ok( "operators/"+opt.label+"."+Obs+".KOP" );
	
	Hk.readLatticeInfo( "operators/"+opt.label+".LAT");	
	Fk.readLatticeInfo( "operators/"+opt.label+".LAT");	
	Ok.readLatticeInfo( "operators/"+opt.label+".LAT");	


	// Get the dimension and distribution of orbitals and supercells
	const int orbPerCell = Hk.Dim();
	const long int numCell = scdim[0]*scdim[1]*scdim[2];

	// The operators are always assumed to be a density function
	// Therefore are rescaled for the equivalent system size
	Ok.Rescale(1.0/kdim[0]/kdim[1]/kdim[2]/numCell);

	// Here we construct the equivalent reciprocal space vector
	// for the particular unit cell we read on K_UnitCell
	// but assuming we expanding into a super cell
	qt::real b[3][3] ;
	for(qt::index i=0;i<3;i++)
	for(qt::index j=0;j<3;j++)
		b[i][j]= Hk.Rec(i,j)/(qt::real)scdim[i];


	// The class LatticeFFTW will be used 
	// to control the parallelization over k
	LatticeFFTW PotFFT;
	PotFFT.Init(MPI_COMM_WORLD,scdim, orbPerCell,opt.label);

	const int 
	loc_dim = PotFFT.LocalDim(), 	//Local dimension dimension of the subspace per process 
	loc_scdim0 = PotFFT.LocalDim0(),//Local supercell dimension from a0 lat-vect
	i0_shift = PotFFT.LocalI0();	//coordinate shift due to the parallelization


	//COMPUTE ESTIMATED MEMORY
	qt::real estimatedMem = 0.0;
	estimatedMem+= 2.0*loc_dim*sizeof(qt::complex);//Increase the estimated memory by loc_dim 
	estimatedMem+= (loc_dim/orbPerCell)*orbPerCell*orbPerCell*sizeof(qt::complex);
	estimatedMem+= 8*loc_dim * sizeof(qt::complex);
	estimatedMem+= MomSize*MomSize*sizeof(qt::complex);
	estimatedMem*=1.5; //scale of memory
	qt::real MaxMem = 1024.*1024.*1024.*20.;

	// Finally we print the options to check everything is ok
	if( world_rank == world_root)
	{
		std::cout<<std::endl;
		std::cout<<"OPERATORS"<<std::endl;
		std::cout<<"The total number of orbitals per cell is: "<<orbPerCell<<std::endl;
		std::cout<<"The Hamiltonian label is "<<Ham<<std::endl;
		std::cout<<"The Force label is " <<For<<std::endl;
		std::cout<<"The Observable label is "<<Obs<<std::endl;
		std::cout<<std::endl<<std::endl;
		std::cout<<"The algorithm is going to compute: "
				 <<MomSize<<"x"<<MomSize<<" Chebyshev moments for the Kubo-Bastin ChebExp formalism"<<std::endl;
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
		std::cout<<" MEMORY_ESTIMATION"<<std::endl;
		std::cout<<"This simulation will use a total ammount of "<< estimatedMem/1024./1024.<<" MB ."<<std::endl;
		std::cout<<"The maximum of allowed memory is  "<< MaxMem/1024./1024.<<" MB ."<<std::endl;
	}
	

	if ( estimatedMem > MaxMem )
	{ 	std::cerr<<"The requested memory is larger that the mem bound"<<std::endl;
		MPI_Finalize(); return 0; }
	// Rescal the operators with their appropiate geometric factors
	// The hamiltonian should be rescale for KPM
	Hk.Rescale(1.0/EnergyScalFact);


	//Before proceding, we should first used the seed for initialize
	//the variables for gsl random vector
	int par_seed = RSeed*(1+ 0.0*getpid())*( world_rank + 2321);    //create a seed that is different for process.
	gsl_rng_env_setup();
	gsl_rng *rng;  // random number generator
	rng = gsl_rng_alloc (gsl_rng_default);     // uses the default
	gsl_rng_set (rng, par_seed);                  // set seed

	//Allocate the moment vector which is of the same size for all processes
	std::vector<qt::complex> mu( MomSize*MomSize,0.0 );


	//Then we read the Potential Operator
	//We will assume that this operator is written also as 
	//a block diagonal operator but in real space,
	sparse::BlockMatrix<qt::complex> Uloc( loc_dim , orbPerCell );
	//Add Electrostatic Disorder
	//Add Structural Disorder
	//Add Adatoms
	for(size_t loc_i0=0; loc_i0 < loc_scdim0 ; loc_i0++ )
	for(size_t i1=0; i1 < scdim[1] 	 ; i1++ )	
	for(size_t i2=0; i2 < scdim[2] 	 ; i2++ )	
	{
		const qt::real W = 2.8;
		const qt::index bidx = ( loc_i0*scdim[1] + i1 )*scdim[2] + i2;
//		for( qt::index o=0; o<orbPerCell; o++ )
//			Uloc(bidx,o,o) = 0.5*W*( 2.0*gsl_rng_uniform(rng) -1.0 )/EnergyScalFact;   
	}
	
	//In addition, we create three BlockDiag Operators 
	//for the hamiltonian, the force, and the operators
	sparse::BlockMatrix<qt::complex>
	Hloc( loc_dim , orbPerCell ),
	Floc( loc_dim , orbPerCell ),
	Oloc( loc_dim , orbPerCell );

	//Allocate the memory for all the necessary arrays (in parallel)
	std::vector<qt::complex> data( 8*loc_dim );
	qt::complex
	*Psi = &data[0*loc_dim],
	*jmL = &data[1*loc_dim],
	*jmR = &data[2*loc_dim],
	*PsiR= &data[3*loc_dim],
	*PsiL= &data[4*loc_dim],
	*jmt = &data[5*loc_dim],
	*jmU = &data[6*loc_dim];
	//Add the size of the vectors

	qt::dimension
	DimBlockMom = (qt::dimension) ((qt::real)(MaxMem-estimatedMem)/loc_dim/sizeof(qt::complex));
	
	if(DimBlockMom > MomSize ) DimBlockMom= MomSize;
	qt::dimension
		BlockNum = ( MomSize+( (qt::integer)DimBlockMom-1) )/DimBlockMom;
	
	if( world_rank == 0)
		std::cout	<<"\nProceding to allocate an block of "
					<<DimBlockMom<<" right Chebyshev vectors"<<std::endl
					<<"In the form of "<<BlockNum<<" blocks"<<std::endl;

	
	std::vector<qt::complex> BlockMomArray( DimBlockMom*loc_dim );
	qt::complex *PsimR[DimBlockMom];
	for( qt::index m=0; m< DimBlockMom ; m++)
		PsimR[m]=&BlockMomArray[m*loc_dim];


	//There are two possible loops, one is the K-point look
	//asociated with the periodic system, and the other is the
	// Nsc-dimensional loop in the momentum space associated with the 
	// the super cell of dimension Nsc. We begin by considering the 
	// Periodic loop, for which we will use the the indexes
	// ik0, ik1, and ik2;
	for(qt::integer ik0=0; ik0 < kdim[0] ; ik0++ )
	for(qt::integer ik1=0; ik1 < kdim[1] ; ik1++ )	
	for(qt::integer ik2=0; ik2 < kdim[2] ; ik2++ )	
	{

		//For each  (ik0,ik1,ik2) Define a k-shift dk in the BZ
		//defined by the periodicity of the unit cell
		// given by:
 		const qt::real
		db[3] = { ik0/(qt::real)kdim[0], ik1/(qt::real)kdim[1]+ ik2/(qt::real)kdim[2] } ;

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
		for(qt::integer loc_i0=0; loc_i0 < loc_scdim0 ; loc_i0++ )
		for(qt::integer i1=0; i1 < scdim[1] 	 ; i1++ )	
		for(qt::integer i2=0; i2 < scdim[2] 	 ; i2++ )	
		{

			const qt::integer bidx = ( loc_i0*scdim[1] + i1 )*scdim[2] + i2;
			const qt::real
			k[3] = {
					(loc_i0+i0_shift + db[0])*b[0][0] + (i1 + db[1])*b[1][0]+ (i2+ db[2])*b[2][0],
					(loc_i0+i0_shift + db[0])*b[0][1] + (i1 + db[1])*b[1][1]+ (i2+ db[2])*b[2][1],
					(loc_i0+i0_shift + db[0])*b[0][2] + (i1 + db[1])*b[1][2]+ (i2+ db[2])*b[2][2]
					};

			Hk.Addk_dependence(k );
			Fk.Addk_dependence(k );
			Ok.Addk_dependence(k );

			Hloc.AddBlock( bidx , Hk.Mat() ) ;
			Floc.AddBlock( bidx , Fk.Mat() ) ;
			Oloc.AddBlock( bidx , Ok.Mat() ) ;
 		}
	
		//Before Starting KPM we define the vector for the trace
		//We make no assumption on anything for the moment
		for( qt::integer r=0; r <  RandNum ; r++ )
		{
			// The random vector will varie from kBZ point
			// to kBZ point to simulate a larger vector
			for( qt::integer  n=0 ; n< loc_dim ; n++)
			{
				qt::real theta_r=  ( 2.0*gsl_rng_uniform(rng) -1.0 )*M_PI;
				Psi[ n ] = qt::complex( cos( theta_r ), sin( theta_r) );
				if( world_rank == 0 && n==10)
					std::cout<<"VECTOR "<<Psi[ n ]<<std::endl;
			}

			//------------------------WE NOW START KPM--------------------/
			Oloc.Multiply(1.0, Psi, 0.0, PsiR);// T_m(H) | O Psi>
			for(qt::integer mBId=0;mBId< BlockNum ; mBId++ )
			{
				for(qt::integer mBR=0;mBR< DimBlockMom ; mBR++ )
				{
					const qt::integer mR= mBId*DimBlockMom+ mBR ;
					// T_m(H) | POW Psi>
					kpm::cheb_evolve( mR , Hloc,PsiR, jmR );
					//We incorporate the effect of the potential
					qt::real bR=2; if ( mR == 0 ) bR=1;
					qt::copy(loc_dim,jmt,PsiR);
					PotFFT.DirectFourierTransform(jmt);
					Uloc.Multiply( bR , jmt , 0 , jmU ); 
					PotFFT.InverseFourierTransform(jmU);
					qt::axpy(loc_dim,fft_norm, jmU,jmR);
					//After evolving add the block to the
					qt::copy(loc_dim,PsimR[mBR],PsiR);
				 }

				//Cheb-evolve the left vector
				qt::copy(loc_dim,PsiL,Psi);
				for(qt::integer mL=0; mL< MomSize ; mL++ )
				{	
					kpm::cheb_evolve( mL , Hloc, PsiL, jmL ); //<Psi T_m(H) |
					//We incorporate the effect of the potential
					qt::real bL=2; if ( mL == 0 ) bL=1;
					qt::copy(loc_dim,jmt,PsiL);
					PotFFT.DirectFourierTransform(jmt);
					Uloc.Multiply( bL , jmt , 0 , jmU ); 
					PotFFT.InverseFourierTransform(jmU);
					qt::axpy(loc_dim,fft_norm, jmU,jmL);

					//multiplu by the force
					Floc.Multiply(1.0, PsiL, 0.0, jmt); //<Psi T_m(H) Op |

					//compute the forces for blocks
					for(qt::integer mBR=0;mBR< DimBlockMom ; mBR++ )
					{
						const qt::integer mR= mBId*DimBlockMom+ mBR ;
						qt::complex mutmp;
						qt::dot(loc_dim, jmt , PsimR[mBR] ,mutmp); //<Psi T_m(H) Op | T_m(H) Pow Psi >
						mu[ mR *MomSize + mL ] += mutmp ;
					}
				}
			}
		}
	}

	//After we finished with all the moments, each cell has parts
	//of the final answer. We will collect all of these in cell 0
	std::vector<qt::complex> mu_Final(MomSize*MomSize,0.0);	
	MPI_Reduce( 
				&mu[0]		,  &mu_Final[0], 
				MomSize*MomSize, MPI_DOUBLE_COMPLEX, 
				MPI_SUM, world_root, MPI_COMM_WORLD
				); 

	// Finally save all the moments
	if( world_rank == world_root )
	{
		std::ofstream mom_file(( opt.label+opt.suffix+".KPMmom").c_str());

		mom_file	<<opt.label+opt.suffix<<" "<<MomSize<<" "<<MomSize <<" " //Label  NumChebMoms
					<<opt.energyShift- opt.energyScale*0.5/opt.cutOff<<" " 	// Lower end of the band
					<<opt.energyShift+ opt.energyScale*0.5/opt.cutOff<<" " 	//Higher end of the band
					<<opt.cutOff<<" "<<std::endl; //Options to be depreciated


		const qt::real hbar= 0.6582119E-15;
		for( qt::index m0=0;m0< MomSize ; m0++ )
		for( qt::index m1=0;m1< MomSize ; m1++ )
		{
			qt::complex muelem= hbar*mu_Final[ m0*MomSize + m1 ];
			mom_file<<m0<<" "<<m1<<" "<<muelem.real()<<" "<<muelem.imag()<<std::endl;
		}
	}


MPI_Finalize();

return 0;}



