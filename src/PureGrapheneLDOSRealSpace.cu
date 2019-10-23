	#include "utilidades.h"
//#include "graphene_cusp.h"
#include "lattice_cusp.h"
#include <cusp/transpose.h>
#include <iostream>
#include <fstream>
#include "kpm_cusp.h"

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
typedef cus0p::coo_matrix<Indice, Scalar, cusp::host_memory>   HCOO;
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
    int Nx, Ny,W,L,D,DD,M, ngpu,nflux,B;	//variables de conteo
	FloatType   Emin, Emax,E0min, E0max,E0;
	FloatType tU,U,p;
    //srand(1121322113);
    srand(time(0));
    
    /************Definimos los parametros iniciales **************************/
	if(argc>=12+1){
		Nx=atoi(argv[1]);                                        //Tamaño del sistema
		Ny=atoi(argv[2]);                                        //Tamaño del sistema
		M=atoi(argv[3]);
        U= atof(argv[4]);
        tU=atof(argv[5]);
		nflux  =atoi(argv[6]);
		p=atof(argv[7]);                                     //Densidad de impurezas
        L=atoi(argv[8]);
        W=atoi(argv[9]);
		E0=atof(argv[10]);                                     //Densidad de impurezas
        //      MachineName=argv[11];
        ngpu=atoi(argv[12]);
    }else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, U, tU,IndexFlux ,p ,L, W ,E0 ,MachineName,GPU"<<std::endl;
        return 0;}
	double tempB=52608*nflux/Nx;
	B=ceil(tempB);
	//std::string String = static_cast<ostringstream*>( &(ostringstream() << B) )->str();
	std::ostringstream ss;ss << B;
	std::cout<<"The magnetic Filed B="<<B<<std::endl;

    /***********Declaracion de directorios de datos de salida **************/
	std::string dospath("data/GrapheneLDOS");
	dospath+="Nx";		dospath+=argv[1];
	dospath+="Ny";		dospath+=argv[2];
	dospath+="M";		dospath+=argv[3];
	dospath+="U";		dospath+=argv[4];
	dospath+="tU";		dospath+=argv[5];
	dospath+="B"     ; dospath+=ss.str();
	dospath+="p";		dospath+=argv[7];
	dospath+="L";		dospath+=argv[8];
	dospath+="W";		dospath+=argv[9];
	dospath+="E0";		dospath+=argv[10];
	dospath+=argv[11] ;
	dospath+=".dat"  ;

	std::string sigmaxxpath("data/GrapheneIMPU");
	sigmaxxpath+="Nx"     ; sigmaxxpath+=argv[1];
	sigmaxxpath+="Ny"     ; sigmaxxpath+=argv[2];
	sigmaxxpath+="M"     ; sigmaxxpath+=argv[3];
	sigmaxxpath+="U"     ; sigmaxxpath+=argv[4];
	sigmaxxpath+="tU"    ; sigmaxxpath+=argv[5];
	sigmaxxpath+="B"     ; sigmaxxpath+=ss.str();
	sigmaxxpath+="p"     ; sigmaxxpath+=argv[7];
	sigmaxxpath+="L"     ; sigmaxxpath+=argv[8];
	sigmaxxpath+="W"     ; sigmaxxpath+=argv[9];
	sigmaxxpath+="E0"	 ; sigmaxxpath+=argv[10];
	sigmaxxpath+=argv[11] ;
	sigmaxxpath+=".dat"  ;

    //Select the gpu device(Default value 0)
   cudaSetDevice(ngpu);
	D 	=	Nx*Ny;
	DD	=	2*D;
    HCOO H (DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins
	graphene::lattice(Nx,Ny,H);
    graphene::TSiteDis(Nx,Ny,H,p,U,tU);
	graphene::Magnetic_Field(Nx,Ny,H,nflux);
	RefineSparse(H);
	//print(H);
	FindingExtrema(H,Emin,Emax,400);
//	FindingExtrema(H,Emin,Emax);
	E0min=0.9*Emin;
	E0max=0.9*Emax;
	LDOS(H,Nx,Ny,M,Emin,Emax,E0,W ,L,dospath);


    return 0;
}
/***************************REEScalando el hamiltoniano******************/
