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

template <typename Matrix, typename FloatType>
        void RAMDOS( Matrix& H,typename Matrix::index_type M,typename Matrix::index_type R, const FloatType Emin,const  FloatType Emax,const FloatType alpha,const std::string outputname, int NE ){
            /********************************Algoritmo de Chebyshev y vectores aleatorios ***********************************/
            typedef typename Matrix::value_type Scalar;
            typedef typename Matrix::index_type Integer;
            typedef typename Matrix::value_type::value_type Floating;
            cusp::array1d	<Scalar,cusp::host_memory> mu_h(M,0);
            cusp::array1d	<Scalar,cusp::device_memory> mu;
			chebyshev::random::delta_moments(H,mu_h,M,R);
			mu=mu_h;
			chebyshev::kernel::jackson_kernel(mu,1.0f/(Floating)R); 
			//Floating lambda=3;chebyshev::kernel::lorentz_kernel(mu,1.0f/(Floating)R,lambda);            
			chebyshev::sum::dosram(mu,Emin,Emax,alpha,(Floating)NE,outputname);   
        }

	template <typename Vector, typename Float>  void dosram(Vector& mu,const Float Emin,const  Float Emax,const Float alpha,const Float NE,const std::string outputname){

		typedef typename Vector::value_type Scalar;
		typedef typename Scalar::value_type Floating;
		Float dEn  =2*alpha/NE;
		thrust::counting_iterator<int> index(1);						//Definimos el inidice del vector como iterador
		Floating a=0;
		cusp::complex<Floating> zero=(Floating)0;	
		Scalar DOS;
		dos_binary<Scalar>		delta_chev(a,mu.size(),a);							//Definimos la operacion binaria del kernel (mu,Index)
		thrust::plus<Scalar>	binary_op1;
		//Calculamos el kernel 
		std::ofstream output_file;   output_file.open(outputname.c_str());

		for(Float En=-alpha;En<=alpha;En=En+dEn){
			a=(2/(((Floating)M_PI)*sqrt(1-En*En)))*(2*alpha/(Emax-Emin));
			delta_chev.En=En;
			delta_chev.a =a;
			DOS=thrust::inner_product(mu.begin(),mu.end(), index,zero,binary_op1,delta_chev);	
			output_file<< 0.5*(En*(Emax-Emin)/alpha+(Emax+Emin))<<" "<<DOS.real()<<std::endl;
			} 	
        output_file.close();
        }


int main(int argc, char *argv[])    			//Para simplificar las corridas elegimos N=argv[1], M=argv[2]
{
    
    /***********Declaracion de apuntadores y variables del programa*********/
    int Nx, Ny,D,DD,M, R;	//variables de conteo
	FloatType Emin, Emax,alpha;
	FloatType U;
	int ngpu;
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
        U= atof(argv[4]);
        R=atoi(argv[5]);
        //      MachineName=argv[6];
        ngpu=atoi(argv[7]);
        SEED=atoi(argv[8]);
    }else{
        std::cout<<"Numero de parametros erroneos, programa finalizado"<<std::endl;
        std::cout<<"Los parametros son:"<<std::endl;
        std::cout<<"Nx, Ny, M, U ,R  ,MachineName,GPU"<<std::endl;
        return 0;}
	srand(SEED);

    /***********Declaracion de directorios de datos de salida **************/
	std::string dataID("");
	dataID+="Nx";		dataID+=argv[1];
	dataID+="Ny";		dataID+=argv[2];
	dataID+="M";		dataID+=argv[3];
	dataID+="U";		dataID+=argv[4];
	dataID+="R";		dataID+=argv[5];
	dataID+=argv[6] ;

	//Creates the file that will hold the Density of states 
	std::string dospath("data/GrapheneDOS");
	dospath.append(dataID);
	dospath+=".dat";

	//Select the gpu device(Default value 0)
	std::cout<<"Setting device "<<ngpu<<std::endl;
   	cudaSetDevice(ngpu);
	D 	=	Nx*Ny;
	DD	=	2*D;
	HCOO H (DD,DD,4*DD); // 10 vecinos 3primeros+6segundos+1onsite. Con 2 spins

	graphene::lattice(Nx,Ny,H);
	graphene::Anderson(Nx,Ny,H,U);
	RefineSparse(H);
	Emin=-3.0;
	Emax= 3.0;
	alpha=0.9;
	chebyshev::Rescale(H,Emin,Emax,alpha);
	cycletime(-1);
	std::cout<<std::endl<<"Calculating DOS"<<std::endl;
	srand(time(0)*t.tv_usec * t.tv_sec * pid);
	chebyshev::random::DOSRAM(H,M,R,Emin,Emax,alpha,dospath,6556+1);
	cycletime(-1);


    return 0;
}
/***************************REEScalando el hamiltoniano******************/
