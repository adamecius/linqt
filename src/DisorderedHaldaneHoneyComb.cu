#include "utilidades.h"
#include "linear_cuda_utilities.h"
//#include "graphene_cusp.h"
#include "lattice_cusp.h"
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
	FloatType   Emin, Emax,E0min, E0max;
	FloatType U,UAB,p,lambda;
    //srand(1121322113);
    srand(time(0));
    
    /************Definimos los parametros iniciales **************************/
	if(argc>=10+1){
		Nx=atoi(argv[1]);                                        //Tamaño del sistema
		Ny=atoi(argv[2]);                                        //Tamaño del sistema
		M=atoi(argv[3]);
        U= atof(argv[4]);
		lambda=atof(argv[5]);
		UAB=atof(argv[6]);
		p=atof(argv[7]);                                     //Densidad de impurezas
        R=atoi(argv[8]);
        //MachineName=argv[9];
        ngpu=atoi(argv[10]);
    }else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, U,lambda, UAB ,p ,R  ,MachineName,GPU"<<std::endl;
        return 0;}

    /***********Declaracion de directorios de datos de salida **************/
	std::string dospath("data/GraphenePolarizedSODOS");
	dospath+="Nx";		dospath+=argv[1];
	dospath+="Ny";		dospath+=argv[2];
	dospath+="M";		dospath+=argv[3];
	dospath+="U";		dospath+=argv[4];
	dospath+="lambda";	dospath+=argv[5];
	dospath+="UAB";	dospath+=argv[6];
	dospath+="p";		dospath+=argv[7];
	dospath+="R";		dospath+=argv[8];
	dospath+=argv[9] ;
	dospath+=".dat"  ;

	std::string sigmaxxpath("data/GraphenePolarizedSOConductivityXX");
	sigmaxxpath+="Nx"     ; sigmaxxpath+=argv[1];
	sigmaxxpath+="Ny"     ; sigmaxxpath+=argv[2];
	sigmaxxpath+="M"     ; sigmaxxpath+=argv[3];
	sigmaxxpath+="U"     ; sigmaxxpath+=argv[4];
	sigmaxxpath+="lambda"; sigmaxxpath+=argv[5];
	sigmaxxpath+="UAB"; sigmaxxpath+=argv[6];
	sigmaxxpath+="p"     ; sigmaxxpath+=argv[7];
	sigmaxxpath+="R"     ; sigmaxxpath+=argv[8];
	sigmaxxpath+=argv[9] ;
	sigmaxxpath+=".dat"  ;

	std::string sigmaxypath("data/GraphenePolarizedSOConductivityXY");
	sigmaxypath+="Nx"     ; sigmaxypath+=argv[1];
	sigmaxypath+="Ny"     ; sigmaxypath+=argv[2];
	sigmaxypath+="M"     ; sigmaxypath+=argv[3];
	sigmaxypath+="U"     ; sigmaxypath+=argv[4];
	sigmaxypath+="lambda"	; sigmaxypath+=argv[5];
	sigmaxypath+="UAB"	; sigmaxypath+=argv[6];
	sigmaxypath+="p"     ; sigmaxypath+=argv[7];
	sigmaxypath+="R"     ; sigmaxypath+=argv[8];
	sigmaxypath+=argv[9] ;
	sigmaxypath+=".dat"  ;

	

    //Select the gpu device(Default value 0)
   cudaSetDevice(ngpu);
	D 	=	Nx*Ny;
	DD	=	2*D;
    	HCOO H (DD,DD,10*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins
	HCOO Vx(DD,DD,10*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins	
	HCOO Vy(DD,DD,10*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins	

	graphene::lattice(Nx,Ny,H);
    graphene::LatticePotential(Nx,Ny,H,(double)1,UAB);
	graphene::Andeson(Nx,Ny,H,U);
	graphene::PolarizedIntrinsecSO(Nx,Ny,H,lambda,p);
	graphene::velocityx(Nx,Ny,H,Vx);
	graphene::velocityy(Nx,Ny,H,Vy);
	RefineSparse(H);
	RefineSparse(Vx);
	RefineSparse(Vy);
//	print(H);
	FindingExtrema(H,Emin,Emax,400);
//	FindingExtrema(H,Emin,Emax);
	E0min=0.97*Emin;
	E0max=0.97*Emax;
//	FindingExtrema(H,Emin,Emax,400);
//	E0min=-3;
//	E0max= 3;	

 	std::cout<<std::endl<<"Calculating DOS"<<std::endl;
	cycletime(-1);
	chebyshev::random::DOS(H,M,R,Emin,Emax,E0min,E0max,dospath,32768); //falta un 8 al final
	cycletime(-1);
	std::cout<<"Calculating SIGMAXY "<<std::endl;
	chebyshev::random::SIGMA(H,Vx,Vy,M,R,Emin,Emax,E0min,E0max,sigmaxypath,16384);
	cycletime(-1);
	std::cout<<"Calculating SIGMAXX "<<std::endl;
	chebyshev::random::SIGMA(H,Vx,Vx,M,R,Emin,Emax,E0min,E0max,sigmaxxpath,16384);//falta un 4 al final


    return 0;
}
/***************************REEScalando el hamiltoniano******************/
	
