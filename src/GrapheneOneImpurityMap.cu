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
    int Nx, Ny, D, DD, M,ngpu;	//variables de conteo
	FloatType Emin, Emax;
	FloatType U,V0,eps0;
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
		M =atoi(argv[3]);
        U = atof(argv[4]);
        V0= atof(argv[5]);
        //      MachineName=argv[6];
        ngpu=atoi(argv[7]);
        SEED=atoi(argv[8]);
		}else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, U, V0, MachineName,GPU, SEED"<<std::endl;
        return 0;}
	srand(SEED);

    /***********Declaracion de directorios de datos de salida **************/
	std::string dataID("");
	dataID+="Nx";		dataID+=argv[1];
	dataID+="Ny";		dataID+=argv[2];
	dataID+="M";		dataID+=argv[3];
	dataID+="U";		dataID+=argv[4];
	dataID+="V0"     ;   dataID+=argv[5];
	dataID+=argv[6] ;

	std::string dosmappath("data/GrapheneHsiteOneImpDOS");
	dosmappath.append(dataID);
	dosmappath+=".dat";
	//Creates the file that will hold the Longitudinal Conductivity
	std::string condmappath("data/GrapheneHsiteOneImpCond");
	condmappath.append(dataID);
	condmappath+=".dat";

	//Select the gpu device(Default value 0)
	std::cout<<"Setting device "<<ngpu<<std::endl;
   	cudaSetDevice(ngpu);


	Scalar V[]={V0,V0,V0,V0,V0,V0}; eps0=1; 
	D 	=	Nx*Ny;
	DD	=	2*D;
	HCOO H (DD,DD,13*DD); // 13 vecinos=1onsite+3primeros+6segundos+3terceros
	HCOO Vx(DD,DD,13*DD); // 13 vecinos=1onsite+3primeros+6segundos+3terceros	
	graphene::lattice(Nx,Ny,H);
	graphene::impurities::Htype(Nx,Ny,H,(int)(Nx/2.0),(int)(Ny/2.0),U,V,eps0);
	graphene::velocityx(Nx,Ny,H,Vx);
	Emin=-4.5;
	Emax= 4.5;
	RefineSparse(H);
	RefineSparse(Vx);
//	linalg::SpectralBounds(H,Emin,Emax,0.0001);
	//print(H);
	chebyshev::Rescale(H,Emin,Emax,1.0);
//	chebyshev::random::DOS(H,M,10,Emin,Emax,(double)0.9,dosmappath,1000);
	chebyshev::realspace::LinearDosConductivity(Nx,Ny,H,Vx,M,Emin,Emax,dosmappath,condmappath,(FloatType)0.01,10.0,10.0);


    return 0;
}
/***************************REEScalando el hamiltoniano******************/
