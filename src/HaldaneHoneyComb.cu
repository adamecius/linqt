#include "utilidades.h"
#include "linear_cuda_utilities.h"
//#include "graphene_cusp.h"
#include "lattice_cusp.h"
#include <iostream>
#include <fstream>
#include "kpm_cusp.h"
#include <sys/time.h>
#include <unistd.h>
#include <sstream>

#define FCOMPLEX
//#define DCOMPLEX
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
    int Nx, Ny,D,DD,M, R,ngpu,En,seed;	//variables de conteo
	FloatType   Emin, Emax,E0min, E0max;
	FloatType U,UAB,p,lambda,alpha;
    int pid=getpid(); // get it as per your OS
	timeval t;
	gettimeofday(&t, NULL);
	std::stringstream ssPID;		//create a stringstream
	ssPID << pid;//add number to the stream
    /************Definimos los parametros iniciales **************************/
	if(argc>=13+1){
		Nx=atoi(argv[1]);                                        //Tamaño del sistema
		Ny=atoi(argv[2]);                                        //Tamaño del sistema
		M=atoi(argv[3]);
        U= atof(argv[4]);
		lambda=atof(argv[5]);
		UAB=atof(argv[6]);
		p=atof(argv[7]);                                     //Densidad de impurezas
        R=atoi(argv[8]);
		alpha=atof(argv[9]);                                     //Densidad de impurezas
        En=atoi(argv[10]);
        //MachineName=argv[11];
        ngpu=atoi(argv[12]);
	seed=atoi(argv[13]);
		}else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, U,lambda, UAB ,p ,R,alpha,En  ,MachineName,GPU, seed"<<std::endl;
        return 0;}
	srand(seed);

    /***********Declaracion de directorios de datos de salida **************/
	std::string dataID("");
	dataID+="Nx";		dataID+=argv[1];
	dataID+="Ny";		dataID+=argv[2];
	dataID+="M";		dataID+=argv[3];
	dataID+="U";		dataID+=argv[4];
	dataID+="lambda";	dataID+=argv[5];
	dataID+="UAB";		dataID+=argv[6];
	dataID+="p";		dataID+=argv[7];
	dataID+="R";		dataID+=argv[8];
	dataID+="alpha";		dataID+=argv[9];
	dataID+="En";		dataID+=argv[10];
		dataID+=argv[11] ;


	//Creates the file that will hold oll the stdout information 
	//Creates the file that will hold the Density of states 
	std::string dospath("data/GrapheneHaldaneDOS");
	dospath.append(dataID);
	dospath+=".dat";
	//Creates the file that will hold the Longitudinal Conductivity
	std::string sigmaxxpath("data/GrapheneHaldaneConductivityXX");
	sigmaxxpath.append(dataID);
	sigmaxxpath+=".dat";
	//Creates the file that will hold the Transverse Conductivity
	std::string sigmaxypath("data/GrapheneHaldaneConductivityXY");
	sigmaxypath.append(dataID);
	sigmaxypath+=".dat";

	

    //Select the gpu device(Default value 0)
   cudaSetDevice(ngpu);
	D 	=	Nx*Ny;
	DD	=	2*D;
    	HCOO H (DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins
	HCOO Vx(DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins	
	HCOO Vy(DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins	

	graphene::lattice(Nx,Ny,H);
//    graphene::LatticePotential(Nx,Ny,H,(double)1,UAB);
//	graphene::TSiteDis(Nx,Ny,H,p,U,(FloatType)0);
//	graphene::Anderson(Nx,Ny,H,U);
//	graphene::HaldaneHoneyComb(Nx,Ny,H,lambda);
	graphene::velocityx(Nx,Ny,H,Vx);
	graphene::velocityy(Nx,Ny,H,Vy);
	RefineSparse(H);
	RefineSparse(Vx);
	RefineSparse(Vy);
	E0min=-En*3.45;
	E0max= En*3.45;
	Emin= E0min;
	Emax= E0max;
	chebyshev::Rescale(H,Emin,Emax,alpha);

	srand(time(0)*t.tv_usec * t.tv_sec * pid);
 	std::cout<<std::endl<<"Calculating DOS"<<std::endl;
	cycletime(-1);
//	chebyshev::random::DOS(H,M,R,Emin,Emax,alpha,dospath,	16384); //falta un 8 al final
	cycletime(-1);
	std::cout<<"Calculating SIGMAXX and SIGMAXY BASTIN"<<std::endl;
	chebyshev::random::SIGMA(H,Vx,Vy,M,R,Emin,Emax,alpha,sigmaxxpath,sigmaxypath,65536+1);
	cycletime(-1);

    return 0;
}
/***************************REEScalando el hamiltoniano******************/
	
