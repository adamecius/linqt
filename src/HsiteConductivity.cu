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
	FloatType U,V0,p,eps0;
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
		M =atoi(argv[3]);
        U = atof(argv[4]);
        V0= atof(argv[5]);
        p = atof(argv[6]);
        R = atof(argv[7]);
        //      MachineName=argv[8];
        ngpu=atoi(argv[9]);
        SEED=atoi(argv[10]);
		}else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, U, V0, p, R, MachineName,GPU, SEED"<<std::endl;
        return 0;}
	srand(SEED);

    /***********Declaracion de directorios de datos de salida **************/
	std::string dataID("");
	dataID+="Nx";		dataID+=argv[1];
	dataID+="Ny";		dataID+=argv[2];
	dataID+="M";		dataID+=argv[3];
	dataID+="U";		dataID+=argv[4];
	dataID+="V0"     ;   dataID+=argv[5];
	dataID+="p"     ;   dataID+=argv[6];
	dataID+="R";		dataID+=argv[7];
	dataID+=argv[8] ;

	//Creates the file that will hold the Density of states 
	std::string dospath("data/GrapheneHsiteImpDOS");
	dospath.append(dataID);
	dospath+=".dat";
	//Creates the file that will hold the Longitudinal Conductivity
	std::string sigmaxxpath("data/GrapheneHsiteImpConductivityXX");
	sigmaxxpath.append(dataID);
	sigmaxxpath+=".dat";
	//Creates the file that will hold the Transverse Conductivity
	std::string sigmaxypath("data/GrapheneHsiteImpConductivityXY");
	sigmaxypath.append(dataID);
	sigmaxypath+=".dat";


	//Select the gpu device(Default value 0)
	std::cout<<"Setting device "<<ngpu<<std::endl;
   	cudaSetDevice(ngpu);
/*As shown in the following draw
*	    V1=(i,j+1)_B ____ V2=(i,j+1)_A
*			/    \
*	    V0=(i,j)_A /      \V3=(i+1,j+1)_B
*		       \      /
*	    V5=(i+1,j)_B\____/V4=(i+1,j)_A
*
* the vector V[6]=(V0,V1,V2,V3,V4,V5) has as components the hybridization amplitude
* of the atoms in the ring within the Anderson Impurity model
*/ 
	//Stype
	Scalar V[]={V0,V0,V0,V0,V0,V0}; eps0=1; Emin=-4.0;Emax= 4.0;
	//Ftype
//	Scalar V[]={V0,-V0,V0,-V0,V0,-V0}; Emin=-4.6;Emax= 4.6; eps0=-1;
	//Ptype
//	Scalar V[]={V0,0,0,-V0,0,0};
	//D1type
//	Scalar V[]={V0,-0.5*V0,-0.5*V0,V0,-0.5*V0,-0.5*V0};
	//D2type
//	Scalar V[]={0,V0,-V0,0,V0,-V0};

	D 	=	Nx*Ny;
	DD	=	2*D;
	HCOO H (DD,DD,13*DD); // 13 vecinos=1onsite+3primeros+6segundos+3terceros
	HCOO Vx(DD,DD,13*DD); // 13 vecinos=1onsite+3primeros+6segundos+3terceros	
	HCOO Vy(DD,DD,13*DD); // 13 vecinos=1onsite+3primeros+6segundos+3terceros
	graphene::lattice(Nx,Ny,H);
	graphene::impurities::HtypeDistribution(Nx,Ny,H,p,(double)0,V,eps0);
//	graphene::Anderson(Nx,Ny,H,U);
	graphene::velocityx(Nx,Ny,H,Vx);
	graphene::velocityy(Nx,Ny,H,Vy);
	linalg::SpectralBounds(H,Emin,Emax,0.001);
	RefineSparse(H);
	RefineSparse(Vx);
	RefineSparse(Vy);
	alpha=0.9;
	chebyshev::Rescale(H,Emin,Emax,alpha);
	std::cout<<std::endl<<"Calculating DOS"<<std::endl;
	srand(time(0)*t.tv_usec * t.tv_sec * pid);
	cycletime(-1);
	chebyshev::random::DOS(H,M,R,Emin,Emax,alpha,dospath,65536+1);
	cycletime(-1);
	std::cout<<"Calculating SIGMAXX and SIGMAXY BASTIN"<<std::endl;
	chebyshev::random::SIGMA(H,Vx,Vy,M,R,Emin,Emax,alpha,sigmaxxpath,sigmaxypath,65536+1);


    return 0;
}
/***************************REEScalando el hamiltoniano******************/
