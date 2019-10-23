#include "utilidades.h"
#include "linear_cuda_utilities.h"
//#include "graphene_cusp.h"
#include "lattice_cusp.h"
#include <cusp/transpose.h>
#include <iostream>
#include <fstream>
#include "kpm_cusp_updated.h"

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
typedef float    Scalar	HCOO Vx(DDD,DDD,S*4*DDD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins	
;
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
	FloatType	alpha,W, Emin, Emax,E0min, E0max;
    //srand(1121322113);
    srand(time(0));
    
    /************Definimos los parametros iniciales **************************/
	if(argc>=7+1){
		Nx=atoi(argv[1]);                                        //Tamaño del sistema
		Ny=atoi(argv[2]);                                        //Tamaño del sistema
		M=atoi(argv[3]);
		W=atof(argv[4]);                                     //Densidad de impurezas
        R=atoi(argv[5]);
		//MachineName=argv[6];
        ngpu=atoi(argv[7]);
    }else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, W ,R  ,MachineName,GPU"<<std::endl;
        return 0;}

    /***********Declaracion de directorios de datos de salida **************/
	std::string dospath("data/GrapheneDOSKBF");
	dospath+="Nx";		dospath+=argv[1];
	dospath+="Ny";		dospath+=argv[2];
	dospath+="M";		dospath+=argv[3];
	dospath+="W";		dospath+=argv[4];
	dospath+="R";		dospath+=argv[5];
	dospath+=argv[6] ;
	dospath+=".dat"  ;

	std::string sigmaxxpath("data/GrapheneConductivityKBF");
	sigmaxxpath+="Nx"     ; sigmaxxpath+=argv[1];
	sigmaxxpath+="Ny"     ; sigmaxxpath+=argv[2];
	sigmaxxpath+="M"     ; sigmaxxpath+=argv[3];
	sigmaxxpath+="W"     ; sigmaxxpath+=argv[4];
	sigmaxxpath+="R"     ; sigmaxxpath+=argv[5];
	sigmaxxpath+=argv[6] ;
	sigmaxxpath+=".dat"  ;

	std::string dospath0("data/GrapheneDOSLIN");
	dospath0+="Nx";		dospath0+=argv[1];
	dospath0+="Ny";		dospath0+=argv[2];
	dospath0+="M";		dospath0+=argv[3];
	dospath0+="W";		dospath0+=argv[4];
	dospath0+="R";		dospath0+=argv[5];
	dospath0+=argv[6] ;
	dospath0+=".dat"  ;

	std::string sigmaxxpath0("data/GrapheneConductivityLIN");
	sigmaxxpath0+="Nx"     ; sigmaxxpath0+=argv[1];
	sigmaxxpath0+="Ny"     ; sigmaxxpath0+=argv[2];
	sigmaxxpath0+="M"     ; sigmaxxpath0+=argv[3];
	sigmaxxpath0+="W"     ; sigmaxxpath0+=argv[4];
	sigmaxxpath0+="R"     ; sigmaxxpath0+=argv[5];
	sigmaxxpath0+=argv[6] ;
	sigmaxxpath0+=".dat"  ;


    //Select the gpu device(Default value 0)
   cudaSetDevice(ngpu);
	D 	=	Nx*Ny;
	DD	=	2*D;
	HCOO H (DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins
	graphene::lattice(Nx,Ny,H);
	graphene::Anderson(Nx,Ny,H,W);
	HCOO V (DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins	
	graphene::velocityx(Nx,Ny,H,V);
	RefineSparse(H);
	RefineSparse(V);
	Emin=-3.3;
	Emax= 3.3;
	alpha=0.9;
	chebyshev::Rescale(H,Emin,Emax,alpha);
	E0min=-1;
	E0max= 1;

	chebyshev::random::DOS(H,M,R,Emin,Emax,alpha,dospath,3111);
	cycletime(-1);
	std::cout<<"Calculating SIGMAXx BASTIN"<<std::endl;
	chebyshev::random::SIGMA(H,V,M,R,Emin,Emax,alpha,sigmaxxpath,3511);
	cycletime(-1);
	chebyshev::random::DosConductivity(H,V,M,R,Emin,Emax,alpha,dospath0,sigmaxxpath0,E0min,E0max,100);

    return 0;
}
/***************************REEScalando el hamiltoniano******************/
