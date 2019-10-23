#include "utilidades.h"
#include "linear_cuda_utilities.h"
//#include "graphene_cusp.h"
#include "lattice_cusp.h"
#include <cusp/transpose.h>
#include <iostream>
#include <fstream>
#include "kpm_cusp.h"
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
	FloatType W;
	int SEED;
    int pid=getpid(); // get it as per your OS
	timeval t;
	gettimeofday(&t, NULL);
	std::stringstream ssPID;		//create a stringstream
	ssPID << pid;//add number to the stream
    
    /************Definimos los parametros iniciales **************************/
	if(argc>=8+1){
		Nx=atoi(argv[1]);                                        //Tamaño del sistema
		Ny=atoi(argv[2]);                                        //Tamaño del sistema
		M=atoi(argv[3]);
		W=atof(argv[4]);
		R=atoi(argv[5]);
        //      MachineName=argv[6];
        ngpu=atoi(argv[7]);
        SEED=atoi(argv[8]);
    }else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, W,R  ,MachineName,GPU"<<std::endl;
        return 0;}
	srand(SEED);

    /***********Declaracion de directorios de datos de salida **************/
	std::string dataID("");
	dataID+="Nx";		dataID+=argv[1];
	dataID+="Ny";		dataID+=argv[2];
	dataID+="M";		dataID+=argv[3];
	dataID+="W";		dataID+=argv[4];
	dataID+="R";		dataID+=argv[5];
	dataID+=argv[6] ;

	//Creates the file that will hold the Density of states 
	std::string dospath("data/GrapheneAndersonDOS");
	dospath.append(dataID);
	dospath+=".dat";
	//Creates the file that will hold the Longitudinal Conductivity
	std::string sigmaxxpath("data/GrapheneAndersonConductivityXX");
	sigmaxxpath.append(dataID);
	sigmaxxpath+=".dat";
	//Creates the file that will hold the Transverse Conductivity
	std::string sigmaxypath("data/GrapheneAndersonConductivityXY");
	sigmaxypath.append(dataID);
	sigmaxypath+=".dat";


	//Select the gpu device(Default value 0)
	std::cout<<"Setting device "<<ngpu<<std::endl;
   	cudaSetDevice(ngpu);
	D 	=	Nx*Ny;
	DD	=	2*D;
	HCOO H (DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins
	HCOO Vx(DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins	
	HCOO Vy(DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins	

	graphene::lattice(Nx,Ny,H);
	graphene::Anderson(Nx,Ny,H,W);
	graphene::velocityx(Nx,Ny,H,Vx);
	graphene::velocityy(Nx,Ny,H,Vy);
	RefineSparse(H);
	RefineSparse(Vx);
	RefineSparse(Vy);
	Emin=-3.05-0.5*W;
	Emax= 3.05+0.5*W;
	alpha=0.1;
	chebyshev::Rescale(H,Emin,Emax,alpha);
	cycletime(-1);
	std::cout<<std::endl<<"Calculating DOS"<<std::endl;
	srand(time(0)*t.tv_usec * t.tv_sec * pid);
	chebyshev::random::DOS(H,M,R,Emin,Emax,alpha,dospath,65536+1);
	cycletime(-1);
	std::cout<<"Calculating SIGMAXX and SIGMAXY BASTIN"<<std::endl;
//	chebyshev::random::SIGMA(H,Vx,Vy,M,R,Emin,Emax,alpha,sigmaxxpath,sigmaxypath,65536+1);


    return 0;
}
/***************************REEScalando el hamiltoniano******************/
