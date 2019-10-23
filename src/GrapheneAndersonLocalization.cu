#include "utilidades.h"
#include "linear_cuda_utilities.h"
#include "graphene_lattice_cusp.h"
#include <cusp/transpose.h>
#include <iostream>
#include <fstream>
#include "kpm_cusp_updated.h"
#include <sys/time.h>
#include <unistd.h>
#include <sstream>
//#define FCOMPLEX
#define DCOMPLEX
//#define FLOAT
//#define DOUBLE

#ifdef FCOMPLEX
cusp::complex<float> I(0.0f,1.0f);
cusp::complex<float> zero(0.0f,0.0f);
typedef float    FloatType;
typedef cusp::complex<float>    Scalar;
typedef int         Indice;
typedef cusp::coo_matrix<Indice, Scalar, cusp::device_memory> DCOO;
typedef cusp::coo_matrix<Indice, Scalar, cusp::host_memory>   HCOO;
#endif

#ifdef DCOMPLEX
cusp::complex<double> I(0.0f,1.0f);
cusp::complex<double> zero(0.0f,0.0f);
typedef double    FloatType;
typedef cusp::complex<double>    Scalar;
typedef int         Indice;
typedef cusp::coo_matrix<Indice, Scalar, cusp::device_memory> DCOO;
typedef cusp::coo_matrix<Indice, Scalar, cusp::host_memory>   HCOO;
#endif


#ifdef FLOAT
typedef float    FloatType;
typedef int      Indice;
typedef cusp::coo_matrix<Indice, Scalar, cusp::device_memory> DCOO;
typedef cusp::coo_matrix<Indice, Scalar, cusp::host_memory>   HCOO;
cusp::complex<float> I(0.0f,1.0f);
float zero=0;
#endif

#ifdef DOUBLE
typedef double    Scalar;
typedef double    FloatType;
typedef int       Indice;
typedef cusp::coo_matrix<Indice, Scalar, cusp::device_memory> DCOO;
typedef cusp::coo_matrix<Indice, Scalar, cusp::host_memory>   HCOO;
cusp::complex<double> I(0.0,1.0);
double zero=0;
#endif


int main(int argc, char *argv[])    			//Para simplificar las corridas elegimos N=argv[1], M=argv[2]
{
    
    /***********Declaracion de apuntadores y variables del programa*********/
    int Nx, Ny,D,DD,M, R,ngpu;	//variables de conteo
	FloatType Emin, Emax,alpha;
	FloatType W,Lx,Ly;
	int SEED;
    int pid=getpid(); // get it as per your OS
	timeval t;
	gettimeofday(&t, NULL);
	std::stringstream ssPID;		//create a stringstream
	ssPID << pid;//add number to the stream
    
    /************Definimos los parametros iniciales **************************/
	if(argc>=10+1){
		Nx=atoi(argv[1]);                                        //Tamaño del sistema
		Ny=atoi(argv[2]);                                        //Tamaño del sistema
		M=atoi(argv[3]);
        W= atof(argv[4]);
        R=atoi(argv[5]);
		Lx=atof(argv[6]);
		Ly=atof(argv[7]);
        //      MachineName=argv[8];
        ngpu=atoi(argv[9]);
        SEED=atoi(argv[10]);
    }else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, W, R,Lx, Ly  ,MachineName,GPU"<<std::endl;
        return 0;}
	srand(SEED);

    /***********Declaracion de directorios de datos de salida **************/
	std::string dataID("");
	dataID+="Nx";		dataID+=argv[1];
	dataID+="Ny";		dataID+=argv[2];
	dataID+="M";		dataID+=argv[3];
	dataID+="W";		dataID+=argv[4];
	dataID+="R";		dataID+=argv[5];
	dataID+="Lx";		dataID+=argv[6];
	dataID+="Ly";		dataID+=argv[7];
	dataID+=argv[8] ;

	//Creates the file that will hold the Density of states 
	std::string dospath("data/GrapheneDOS");
	dospath.append(dataID);
	dospath+=".dat";

	//Creates the file that will hold the Density of states 
	std::string olddospath("data/GrapheneOLDDOS");
	olddospath.append(dataID);
	olddospath+=".dat";


	//Creates the file that will hold the Longitudinal Conductivity
	std::string spec_sigmaxx_path("data/GrapheneSpectralConductivityXX");
	spec_sigmaxx_path.append(dataID);
	spec_sigmaxx_path+=".dat";

	//Creates the file that will hold the Longitudinal Conductivity
	std::string spec_sigmaxy_path("data/GrapheneSpectralConductivityXY");
	spec_sigmaxy_path.append(dataID);
	spec_sigmaxy_path+=".dat";


	//Creates the file that will hold the Longitudinal Conductivity
	std::string sigmaxxpath("data/GrapheneConductivityXX");
	sigmaxxpath.append(dataID);
	sigmaxxpath+=".dat";

	//Creates the file that will hold the Density of states 
	std::string dosmappath("data/GrapheneDOSMAP");
	dosmappath.append(dataID);
	dosmappath+=".dat";
	//Creates the file that will hold the Longitudinal Conductivity
	std::string sigmaxxmappath("data/GrapheneConductivityXXMAP");
	sigmaxxmappath.append(dataID);
	sigmaxxmappath+=".dat";

	//Creates the file that will hold the Transverse Conductivity
	std::string sigmaxypath("data/GrapheneMagneticConductivityXY");
	sigmaxypath.append(dataID);
	sigmaxypath+=".dat";


	//Select the gpu device(Default value 0)
	std::cout<<"Setting device "<<ngpu<<std::endl;
   	cudaSetDevice(ngpu);
	D 	=	Nx*Ny;
	DD	=	2*D;
	HCOO H (DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins
	HCOO Vx(DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins	

	graphene::lattice(Nx,Ny,H);
	graphene::Anderson(Nx,Ny,H,W);
	graphene::velocityx(Nx,Ny,H,Vx);
	RefineSparse(H);
	RefineSparse(Vx);
	Emin=-4.15;
	Emax= 4.15;
	alpha=0.9;
	chebyshev::Rescale(H,Emin,Emax,alpha);
	cycletime(-1);
	std::cout<<"Calculating OLD DOS"<<std::endl;
    chebyshev::random::DOS(H,M,R,Emin,Emax,alpha,	olddospath,1001);
	cycletime(-1);
	cycletime(-1);
	std::cout<<"Calculating Spectral Conductivity"<<std::endl;
	//chebyshev::random::SIGMA(H,Vx,Vx,M,R,Emin,Emax,alpha,spec_sigmaxx_path,spec_sigmaxy_path,65341);
	cycletime(-1);
	cycletime(-1);
	std::cout<<"Calculating DOS and ZERO TEMP Conductivity"<<std::endl;
	//chebyshev::random::LinearDosConductivity(H,Vx,M,R,Emin,Emax,alpha,dospath,sigmaxxpath,-3.0,3.0,100);
	cycletime(-1);
	cycletime(-1);
	std::cout<<"Calculating DOS and ZERO TEMP Conductivity maps"<<std::endl;
	chebyshev::realspace::LinearDosConductivity(Nx,Ny,H,Vx,M,Emin,Emax,alpha,dosmappath,sigmaxxmappath,(FloatType)0.01,sqrt(3)*(Lx+0.5*Ly),1.5*Ly);
	cycletime(-1);


    return 0;
}
/***************************REEScalando el hamiltoniano******************/
