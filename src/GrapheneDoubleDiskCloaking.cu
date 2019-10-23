#include "utilidades.h"
#include "linear_cuda_utilities.h"
//#include "graphene_cusp.h"
#include "graphene_lattice_cusp.h"
//#include "lattice_cusp.h"
#include <cusp/transpose.h>
#include <iostream>
#include <fstream>
#include "kpm_cusp_updated.h"
//#include "kpm_cusp.h"

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
    int Nx, Ny, D, DD, M, R,ngpu;	//variables de conteo
	FloatType Emin, Emax,alpha;
	FloatType W, Umin, Umax, Rmin, Rmax;
	int SEED;
    int pid=getpid(); // get it as per your OS
	timeval t;
	gettimeofday(&t, NULL);
	std::stringstream ssPID;		//create a stringstream
	ssPID << pid;//add number to the stream
    
    /************Definimos los parametros iniciales **************************/
	if(argc>=12+1){
		Nx	=atoi(argv[1]);                                        //Tamaño del sistema
		Ny	=atoi(argv[2]);                                        //Tamaño del sistema
		M 	=atoi(argv[3]);
        W 	=atof(argv[4]);
        Umin=atof(argv[5]);
        Rmin=atof(argv[6]);
        Umax=atof(argv[7]);
        Rmax=atof(argv[8]);
        R = atof(argv[9]);
        //      MachineName=argv[10];
        ngpu=atoi(argv[11]);
        SEED=atoi(argv[12]);
		}else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, W, Umin, Rmin, Umax, Rmax, R, MachineName,GPU, SEED"<<std::endl;
        return 0;}
	srand(SEED);

    /***********Declaracion de directorios de datos de salida **************/
	std::string dataID("");
	dataID+="Nx";		dataID+=argv[1];
	dataID+="Ny";		dataID+=argv[2];
	dataID+="M";		dataID+=argv[3];
	dataID+="W";		dataID+=argv[4];
	dataID+="Umin";		dataID+=argv[5];
	dataID+="Rmin";		dataID+=argv[6];
	dataID+="Umax";		dataID+=argv[7];
	dataID+="Rmax";		dataID+=argv[8];
	dataID+="R";		dataID+=argv[9];
	dataID+=argv[8] ;

	//Creates the file that will hold the Density of states 
	std::string dospath("data/GrapheneDoubleDiskDOS");
	dospath.append(dataID);
	dospath+=".dat";
	//Creates the file that will hold the Longitudinal Conductivity
	std::string sigmaxxpath("data/GrapheneDoubleDiskCondXX");
	sigmaxxpath.append(dataID);
	sigmaxxpath+=".dat";
	//Creates the file that will hold the Density of states 
	std::string dosmappath("data/GrapheneDoubleDiskDOSMAP");
	dosmappath.append(dataID);
	dosmappath+=".dat";
	//Creates the file that will hold the Longitudinal Conductivity
	std::string sigmaxxmappath("data/GrapheneDoubleDiskCondXXMAP");
	sigmaxxmappath.append(dataID);
	sigmaxxmappath+=".dat";
	//Select the gpu device(Default value 0)
	std::cout<<"Setting device "<<ngpu<<std::endl;
   	cudaSetDevice(ngpu);


	D 	=	Nx*Ny;
	DD	=	2*D;
	HCOO H (DD,DD,4*DD); // 13 vecinos=1onsite+3primeros+6segundos+3terceros
	HCOO Vx(DD,DD,4*DD); // 13 vecinos=1onsite+3primeros+6segundos+3terceros	

	graphene::lattice(Nx,Ny,H);
	graphene::Anderson(Nx,Ny,H,W);
	graphene::potentials::CenteredDoubleDisk(Nx,Ny,H,Umin,Rmin,Umax,Rmax);
	graphene::velocityx(Nx,Ny,H,Vx);
	Emin=-3.13; Emax= 3.13;
	linalg::SpectralBounds(H,Emin,Emax,0.001);
	RefineSparse(H);
	RefineSparse(Vx);
	alpha=1.0;
	chebyshev::Rescale(H,Emin,Emax,alpha);
	std::cout<<std::endl<<"Calculating DOS"<<std::endl;
	srand(time(0)*t.tv_usec * t.tv_sec * pid);
	cycletime(-1);
	chebyshev::random::DOS(H,M,R,Emin,Emax,alpha,dospath,1000+1);
	chebyshev::realspace::LinearDosConductivity(Nx,Ny,H,Vx,M,Emin,Emax,dosmappath,sigmaxxmappath,(FloatType)0.01,sqrt(3)*(Nx+0.5*Ny),1.5*Ny);
//	cycletime(-1);
//	std::cout<<"Calculating SIGMAXX and SIGMAXY BASTIN"<<std::endl;
//	chebyshev::random::SIGMA(H,Vx,Vy,M,R,Emin,Emax,alpha,sigmaxxpath,sigmaxypath,65536+1);


    return 0;
}
/***************************REEScalando el hamiltoniano******************/
