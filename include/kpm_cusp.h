#pragma once

#ifndef KPM_CUSP_H
#define KPM_CUSP_H
#endif

#include <thrust/functional.h>
#include "utilidades.h"
#include <thrust/inner_product.h>


namespace chebyshev{


/*******************************Finding Extrema****************************************/
template <typename Matrix,typename FloatType>	void Rescale( Matrix& H,const FloatType& Emin,const FloatType& Emax)
{
    
    typedef typename Matrix::value_type Scalar;
    typedef typename Matrix::index_type Integer;
    typedef typename Matrix::memory_space MemorySpace;
    typedef typename Matrix::value_type::value_type Floating;

    cusp::coo_matrix<Integer, Scalar, cusp::device_memory> GH;
	Integer DD=H.num_cols;					//Dimensions for the matrix DD*DD
	GH = H;
    cusp::coo_matrix<Integer, Scalar, cusp::device_memory> Em(DD,DD,DD);

    std::cout<<"We will rescale the hamiltonian using Emin="<<Emin<<" Emax="<<Emax<<std::endl;
    Floating dE=(Emax-Emin)*0.5;
    for(int i=0;i<DD;i++){
        Em.values[i]        =(Emax+Emin)*0.5;
        Em.column_indices[i]=i;
        Em.row_indices[i]   =i;
    }
	cusp::subtract(GH, Em, GH);
    cusp::VectorTransform::matrixscale(1/dE, GH);
    H=GH;    
}

template <typename Matrix,typename FloatType>	void Rescale( Matrix& H,const FloatType& Emin,const FloatType& Emax,const FloatType& alpha)
{
    
    typedef typename Matrix::value_type Scalar;
    typedef typename Matrix::index_type Integer;
    typedef typename Matrix::memory_space MemorySpace;
    typedef typename Matrix::value_type::value_type Floating;

    cusp::coo_matrix<Integer, Scalar, cusp::device_memory> GH;
	Integer DD=H.num_cols;					//Dimensions for the matrix DD*DD
	GH = H;
    cusp::coo_matrix<Integer, Scalar, cusp::device_memory> Em(DD,DD,DD);

    std::cout<<"We will rescale the hamiltonian using Emin="<<Emin<<" Emax="<<Emax<<std::endl;
    Floating dE=(Emax-Emin)*0.5/alpha;
    for(int i=0;i<DD;i++){
        Em.values[i]        =(Emax+Emin)*0.5;
        Em.column_indices[i]=i;
        Em.row_indices[i]   =i;
    }
	cusp::subtract(GH, Em, GH);
    cusp::VectorTransform::matrixscale(1/dE, GH);
    H=GH;    
}

	
namespace kernel{
	
	template<typename Scalar> struct jackson_binary : public thrust::binary_function<Scalar,Scalar,Scalar>{
            typedef typename Scalar::value_type Floating;
			Floating 	M,a;
			jackson_binary(Floating _M,Floating _a) : M(_M), a(_a){}	
			__host__ __device__ Scalar operator()(const Scalar &x, const Floating &m) const {
			//REMEMBER m is a counter and start in 1
			return a*x*(1.0f/(M+1.0f))*((M-m+(Floating)2)*cusp::cos((Floating)M_PI*(m-(Floating)1)/(M+1.0f))+cusp::sin((Floating)M_PI*(m-(Floating)1)/(M+1.0f))*cusp::cos((Floating)M_PI/(M+1.0f))/cusp::sin((Floating)M_PI/(M+1.0f)));}
		}; // end plus

	template<typename Scalar> struct lorentz_binary : public thrust::binary_function<Scalar,Scalar,Scalar>{
            typedef typename Scalar::value_type Floating;
			Floating 	M,a,lambda;
			lorentz_binary(Floating _M,Floating _a,Floating _lambda) : M(_M), a(_a), lambda(_lambda){}	
			__host__ __device__ Scalar operator()(const Scalar &x, const Floating &m) const {
			//REMEMBER m is a counter and start in 1
			return a*x*sinh(lambda*(1-(m-1)/M))/sinh(lambda);}
		}; // end plus


	template<typename Scalar> struct jackson2D_binary : public thrust::binary_function<Scalar,Scalar,Scalar>{
            typedef typename Scalar::value_type Floating;
			Floating 	M,a;
			jackson2D_binary(Floating _M,Floating _a) : M(_M), a(_a){}	
			__host__ __device__ Scalar operator()(const Scalar &x, const Floating &m) const {
			//REMEMBER m is a counter and start in 1
			return a*x*
			(1.0f/(M+1.0f))*((M-floor((m-1)/M)+1				 )*cusp::cos((Floating)M_PI*floor((m-1)/M)				 /(M+1.0f))+cusp::sin((Floating)M_PI*floor((m-1)/M)				  /(M+1.0f))*cusp::cos((Floating)M_PI/(M+1.0f))/cusp::sin((Floating)M_PI/(M+1.0f)))*
			(1.0f/(M+1.0f))*((M-(Floating)(((int)(m-1))%((int)M))+1)*cusp::cos((Floating)M_PI*(Floating)(((int)(m-1))%((int)M))/(M+1.0f))+cusp::sin((Floating)M_PI*(Floating)(((int)(m-1))%((int)M))/(M+1.0f))*cusp::cos((Floating)M_PI/(M+1.0f))/cusp::sin((Floating)M_PI/(M+1.0f)));}
		}; // end plus


	template <typename Vector,typename Float>  void jackson_kernel(Vector& x, Float a){
		typedef typename Vector::value_type T;
		thrust::counting_iterator<int> index(1);						//Definimos el inidice del vector como iterador
		jackson_binary<T> kernel_bin(x.size(),a);							//Definimos la operacion binaria del kernel (mu,Index)
		//Calculamos el kernel 
		thrust::transform(x.begin(),x.end(),index,x.begin(),kernel_bin); 	
        }


	template <typename Vector,typename Float>  void lorentz_kernel(Vector& x, Float a,Float lambda){
		typedef typename Vector::value_type T;
		thrust::counting_iterator<int> index(1);						//Definimos el inidice del vector como iterador
		lorentz_binary<T> kernel_bin(x.size(),a,lambda);							//Definimos la operacion binaria del kernel (mu,Index)
		//Calculamos el kernel 
		thrust::transform(x.begin(),x.end(),index,x.begin(),kernel_bin); 	
        }


	template <typename Vector,typename Float>  void jackson2D_kernel(Vector& mu,  Float a){
		typedef typename Vector::value_type T;
		thrust::counting_iterator<int> index(1);						//Definimos el inidice del vector como iterador
		jackson2D_binary<T> kernel_bin(sqrt(mu.size()),a);							//Definimos la operacion binaria del kernel (mu,Index)
		//Calculamos el kernel 
		thrust::transform(mu.begin(),mu.end(),index,mu.begin(),kernel_bin); 	
        }

	
	}

namespace sum{
	
	template<typename Scalar> struct dos_binary : public thrust::binary_function<Scalar,Scalar,Scalar>{
            typedef typename Scalar::value_type Floating;
			Floating 	M,a,En;
			dos_binary(Floating _En,Floating _M,Floating _a) : En(_En), M(_M), a(_a){}	
			__host__ __device__ Scalar operator()(const Scalar &mu, const Scalar &m) const {
			return a*mu*cusp::cos((m-(Floating)1)*cusp::acos(En));}
		}; // end plus

	template <typename Vector, typename Float>  void dos(Vector& mu,const Float Emin,const  Float Emax,const Float alpha,const Float NE,const std::string outputname){

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

	template<typename Scalar> struct Spectral_DCConductivity_binary : public thrust::binary_function<Scalar,Scalar,Scalar>{
		typedef typename Scalar::value_type Floating;
		Floating	En;				//member_variable
		int 	M;
		Floating  a;
		Scalar	     I;			//member_variable
		Spectral_DCConductivity_binary(Floating _En,int _M,Floating _a) : En(_En), M(_M), a(_a), I(0,1)  {}	
		__host__ __device__ Scalar operator()(const Scalar &mu, const Floating &m) const {
		//floor((m-1)/(Floating)M)	((Floating)((int)(m-1)%M))
		return mu*a*(
				cos(floor((m-1)/(Floating)M  )*acos(En))*cusp::exp(-I*((Floating)((int)(m-1)%M))*acos(En))*(En+I*((Floating)((int)(m-1)%M))*sqrt(1-En*En))+
				cos(((Floating)((int)(m-1)%M))*acos(En))*cusp::exp(+I*floor((m-1)/(Floating)M  )*acos(En))*(En-I*floor((m-1)/(Floating)M  )*sqrt(1-En*En)));}
		}; // end plus

	template <typename Vector,typename Float >  void Spectral_DCConductivity(Vector& mu,const Float Emin,const  Float Emax,const Float alpha, const Float NE,const std::string outputname){

		typedef typename Vector::value_type Scalar;
		typedef typename Scalar::value_type Floating;
		Floating dEn  =2*alpha/NE;
		thrust::counting_iterator<int> index(1);						//Definimos el inidice del vector como iterador
		Floating a=0;
		Floating hbar=0.65821;
		Floating pi=M_PI;

		int M=sqrt(mu.size());
		cusp::complex<Floating> zero=(Floating)0;	
		Scalar SIGMA;
		Spectral_DCConductivity_binary<Scalar>		delta_chev(a,M,a);							//Definimos la operacion binaria del kernel (mu,Index)
		thrust::plus<Scalar>	binary_op1;
		//Calculamos el kernel 
		std::ofstream output_file;   output_file.open(outputname.c_str());
		for(Floating En=-alpha;En<=alpha;En=En+dEn){
			a=2*hbar*hbar*pi/(pow((Emax-Emin)/(2*alpha),3)*pow(1-En*En,2));
			delta_chev.En=En;
			delta_chev.a =a;
			SIGMA=thrust::inner_product(mu.begin(),mu.end(), index,zero,binary_op1,delta_chev);	
			output_file<< 0.5*(En*(Emax-Emin)/alpha+(Emax+Emin))<<" "<<SIGMA.real()<<std::endl;
			} 	
        output_file.close();
        }
        
	
	}

namespace random{
	
	template <typename Matrix, typename Vector>
        void delta_moments( Matrix& H, Vector& mu,int M, int R){ //This method will calculated the chebyshev moment using random vector method.
            /********************************Algoritmo de Chebyshev y vectores aleatorios ***********************************/
            typedef typename Matrix::value_type Scalar;
            typedef typename Matrix::index_type Integer;
            typedef typename Matrix::value_type::value_type Floating;
            
			cusparseHandle_t	cusparse_handle;		//Handle of cusparse
			cusparseMatDescr_t	cusparse_descr;			//Cusparse's Matrix Descriptor  
			cusparseHybMat_t	GH;						//Cusparse's Hybrid Matrix
			cublasHandle_t		cublas_handle;			//Handle of cublas
			cusparseCreate			(&cusparse_handle);	//Initialization of the cusparse's handler
			cusparseCreateMatDescr	(&cusparse_descr);	//Initialization of the descriptor			
			cusparseCreateHybMat	(&GH);				//Initialization of the Hybrid Matrix
			cublasCreate			(&cublas_handle);	//Initialization of the cublas's handler  
            Integer DD=H.num_cols;					//Dimensions for the matrix DD*DD
			cusparse::convertX2Hyb(H,DD,cusparse_handle,cusparse_descr,GH);
            cusp::array1d	<Scalar     , 		cusp::device_memory> jn0(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> jn1(DD,0);
			cusp::array1d	<Scalar     ,       cusp::host_memory> ranx(DD,0);
			Scalar* pjn0=thrust::raw_pointer_cast(jn0.data());
			Scalar* pjn1=thrust::raw_pointer_cast(jn1.data());
            Scalar alpha=	 2.0f;
            Scalar beta =	-1.0f;
            Scalar zero =	 0.0f;
            Scalar one	=	 1.0f;
					M	=	(int)((Floating)M/(Floating)2);
            Scalar mutemp,mu1;
            Floating munorm;
			//Define the first moment mu0=1/2
            mu[0]=((Floating)0.5)*((Scalar)R);
            for(Integer r=0;r<R;r++){
					//Creates the gaussian unormalized random vector on host
                for(Integer i=0;i<DD;i++) ranx[i]=cusp::sqrt(-((Floating)2)*cusp::log((Scalar)rand()/(Scalar)RAND_MAX))*cusp::cos((((Floating)2)*((Floating)M_PI))*(Scalar)rand()/(Scalar)RAND_MAX);
				//Pass the random vector from the host to the device name it jn0 and normalize it
                jn0=ranx;									
                cublas::VectorTransform::normalize<Floating>(cublas_handle,DD,pjn0);//Here we normalize j0 and then we pass it to jn
				//Define jn1=H|jn0>
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pjn0,&zero,pjn1);
				//Calculate the moment mu1=<jn0|H|jn0>=<jn0|jn1>
				cublas::VectorTransform::dot(cublas_handle,DD,pjn0,pjn1	,&mu1);
				mu[1]=mu[1]+mu1;
				//Calculate the moment mu2=2<jn1|jn1>-1= ||jn1>|^2 -1 
				cublas::VectorTransform::nrm2(cublas_handle,DD,pjn1,&munorm);
				mu[2*1+0]=mu[2*1+0]+2*munorm*munorm-(Floating)1;
                for(Integer m=2;m<M;m++){
       	            cycletime(R*(M-2));
					// Redefine the vector |jn0>=2H|jn1>-|jn0> 
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjn1,&beta,pjn0);
					// Calculate the moments mu3,mu5,mu7,...,mu_(2m-1)  
					// using mu_(2m-1)=2<jn0|jn1>-mu1 
					cublas::VectorTransform::dot(cublas_handle,DD,pjn1,pjn0,&mutemp);
					mu[2*m-1]=mu[2*m-1]+((Floating)2)*mutemp		 -mu1;
					// Calculate the moments mu4,mu6,mu7,...,mu_(2m)  
					// using mu_(2m  )=2<jn1|jn1>-mu0 
					cublas::VectorTransform::nrm2(cublas_handle,DD,pjn0,&munorm);
					mu[2*m+0]=mu[2*m+0]+ 2*munorm*munorm-(Floating)1;
					// Exchange jn0<-->jn1 
					cublas::VectorTransform::swap(cublas_handle,DD,pjn0,pjn1);
       	            cycletime(R*(M-2));
                }
            }
			cusparseDestroy			(cusparse_handle);
			cusparseDestroyMatDescr	(cusparse_descr);
			cusparseDestroyHybMat	(GH);
			cublasDestroy			(cublas_handle);
        }

	template <typename Matrix, typename Vector>
        void GreenFun2D_moments( Matrix& H,Matrix& Vx,Matrix& Vy, Vector& muxx, Vector& muxy,int M, int R){ //This method will calculated the chebyshev moment using random vector method.
            /********************************Algoritmo de Chebyshev y vectores aleatorios ***********************************/

            typedef typename Matrix::value_type Scalar;
            typedef typename Matrix::index_type Integer;
            typedef typename Matrix::value_type::value_type Floating;

			cusparseHandle_t	cusparse_handle;		//Handle of cusparse
			cusparseMatDescr_t	cusparse_descr;			//Cusparse's Matrix Descriptor  
			cusparseHybMat_t	GH;						//Cusparse's Hybrid Matrix
			cusparseHybMat_t	GVx;						//Cusparse's Hybrid Matrix
			cusparseHybMat_t	GVy;						//Cusparse's Hybrid Matrix
			cublasHandle_t		cublas_handle;			//Handle of cublas
			cusparseCreate			(&cusparse_handle);	//Initialization of the cusparse's handler
			cusparseCreateMatDescr	(&cusparse_descr);	//Initialization of the descriptor			
			cusparseCreateHybMat	(&GH);				//Initialization of the Hybrid Matrix
			cusparseCreateHybMat	(&GVx);				//Initialization of the Hybrid Matrix
			cusparseCreateHybMat	(&GVy);				//Initialization of the Hybrid Matrix
			cublasCreate			(&cublas_handle);	//Initialization of the cublas's handler

 			//Definitions of integers variables
			Integer DD=H.num_cols;  			
			// Transformations from host to the device matrices
			cusparse::convertX2Hyb(H ,DD,cusparse_handle,cusparse_descr,GH);
			cusparse::convertX2Hyb(Vx,DD,cusparse_handle,cusparse_descr,GVx);
			cusparse::convertX2Hyb(Vy,DD,cusparse_handle,cusparse_descr,GVy);
            //Definitions of the device and hosts vectors
            cusp::array1d	<Scalar     ,       cusp::host_memory>   	HOSTphi(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	Phitxx(DD,0);   
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	Phitxy(DD,0);   
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	Phi(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jn0(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jn1(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jm0(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jm1(DD,0);
			Scalar* pPhi =thrust::raw_pointer_cast(Phi.data());
			Scalar* pPhitxx=thrust::raw_pointer_cast(Phitxx.data());
			Scalar* pPhitxy=thrust::raw_pointer_cast(Phitxy.data());
			Scalar* pjn0 =thrust::raw_pointer_cast(jn0.data());
			Scalar* pjn1 =thrust::raw_pointer_cast(jn1.data());
			Scalar* pjm0 =thrust::raw_pointer_cast(jm0.data());
			Scalar* pjm1 =thrust::raw_pointer_cast(jm1.data());

            Scalar alpha=	 2.0f;
            Scalar beta =	-1.0f;
            Scalar zero =	 0.0f;
            Scalar one	=	 1.0f;
            Scalar mutemp;
			for(Integer r=0;r<R;r++){
				for(Integer i=0;i<DD;i++) 
			//HOSTphi[i]=0;	
HOSTphi[i]=( (Floating)rand()/(Floating)RAND_MAX  -0.5);			
//HOSTphi[100]=1;
				Phi=HOSTphi;
	            cublas::VectorTransform::normalize<Floating>(cublas_handle,DD,pPhi);//Normalize the vectors
				// m=0     n=0//
				//Ket
				jn0=Phi;
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GVx,pjn0,&zero,pPhitxx);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GVy,pjn0,&zero,pPhitxy);
				//Bra
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GVx,pPhi,&zero,pjm0);
				cycletime((Floating)(M*M*R));
                cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxx,&mutemp);//conjugate first vector
				muxx[0*M+0]=muxx[0*M+0]+ mutemp;
                cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);//conjugate first vector
				muxy[0*M+0]=muxy[0*M+0]+ mutemp;
				cycletime((Floating)(M*M*R));
				// m=1     n=0//
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pjm0,&zero,pjm1);
				cycletime((Floating)(M*M*R));
                cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhitxx,&mutemp);
				muxx[1*M+0]=muxx[1*M+0]+ mutemp;
                cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhitxy,&mutemp);
				muxy[1*M+0]=muxy[1*M+0]+ mutemp;
				cycletime((Floating)(M*M*R));
				// all n's left  and m=0//
				for(int m=2;m<M;m++){
					cycletime((Floating)(M*M*R));
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxx,&mutemp);
					muxx[m*M+0]=muxx[m*M+0]+mutemp;
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);
					muxy[m*M+0]=muxy[m*M+0]+mutemp;
					cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
       				cycletime((Floating)(M*M*R));
					}
				//We now goes to the following m=1, and restart n=0;
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH ,pjn0,&zero,pjn1);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GVx,pjn1,&zero,pPhitxx);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GVy,pjn1,&zero,pPhitxy);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GVx,pPhi,&zero,pjm0);
				cycletime((Floating)(M*M*R));
				cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxx,&mutemp);
				muxx[0*M+1]=muxx[0*M+1]+ mutemp;
				cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);
				muxy[0*M+1]=muxy[0*M+1]+ mutemp;
				cycletime((Floating)(M*M*R));
				// n=1     m=1//
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH ,pjm0,&zero,pjm1);
				cycletime((Floating)(M*M*R));
				cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhitxx,&mutemp);				
				muxx[1*M+1]=muxx[1*M+1]+ mutemp;
				cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhitxy,&mutemp);				
				muxy[1*M+1]=muxy[1*M+1]+ mutemp;
				cycletime((Floating)(M*M*R));
				// all m's left  and n=1//
				for(int m=2;m<M;m++){
					cycletime((Floating)(M*M*R));
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxx,&mutemp);
					muxx[m*M+1]=muxx[m*M+1]+mutemp;
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);
					muxy[m*M+1]=muxy[m*M+1]+mutemp;
					cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
					cycletime((Floating)(M*M*R));
					}
				//We now iterate over all n's and restart m;
				for(int n=2;n<M;n++){
					cycletime((Floating)(M*M*R));
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH ,pjn1,&beta,pjn0);
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GVx,pjn0,&zero,pPhitxx);
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GVy,pjn0,&zero,pPhitxy);
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GVx,pPhi,&zero,pjm0);
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxx,&mutemp);
					muxx[0*M+n]=muxx[0*M+n]+mutemp;
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);
					muxy[0*M+n]=muxy[0*M+n]+mutemp;
					cycletime((Floating)(M*M*R));
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GH ,pjm0,&zero,pjm1);
					cycletime((Floating)(M*M*R));
					cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhitxx,&mutemp);
					muxx[1*M+n]=muxx[1*M+n]+ mutemp;
					cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhitxy,&mutemp);
					muxy[1*M+n]=muxy[1*M+n]+ mutemp;
					cycletime((Floating)(M*M*R));
					for(int m=2;m<M;m++){
						cycletime((Floating)(M*M*R));
						cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
						cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxx,&mutemp);
						muxx[m*M+n]=muxx[m*M+n]+mutemp;
						cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);
						muxy[m*M+n]=muxy[m*M+n]+mutemp;
						cycletime((Floating)(M*M*R));
						cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
					}					
				cublas::VectorTransform::swap(cublas_handle,DD,pjn0,pjn1);
				}
			}

        };

	template <typename Matrix, typename Vector>
        void GreenFun2D_moments( Matrix& H,Matrix& V,Vector& mu,int M, int R){ //This method will calculated the chebyshev moment using random vector method.
            /********************************Algoritmo de Chebyshev y vectores aleatorios ***********************************/

            typedef typename Matrix::value_type Scalar;
            typedef typename Matrix::index_type Integer;
            typedef typename Matrix::value_type::value_type Floating;

			cusparseHandle_t	cusparse_handle;		//Handle of cusparse
			cusparseMatDescr_t	cusparse_descr;			//Cusparse's Matrix Descriptor  
			cusparseHybMat_t	GH;						//Cusparse's Hybrid Matrix
			cusparseHybMat_t	GV;						//Cusparse's Hybrid Matrix
			cublasHandle_t		cublas_handle;			//Handle of cublas
			cusparseCreate			(&cusparse_handle);	//Initialization of the cusparse's handler
			cusparseCreateMatDescr	(&cusparse_descr);	//Initialization of the descriptor			
			cusparseCreateHybMat	(&GH);				//Initialization of the Hybrid Matrix
			cusparseCreateHybMat	(&GV);				//Initialization of the Hybrid Matrix
			cublasCreate			(&cublas_handle);	//Initialization of the cublas's handler
 			//Definitions of integers variables
			Integer DD=H.num_cols;  			
			// Transformations from host to the device matrices
			cusparse::convertX2Hyb(H ,DD,cusparse_handle,cusparse_descr,GH);
			cusparse::convertX2Hyb(V,DD,cusparse_handle,cusparse_descr,GV);
            //Definitions of the device and hosts vectors
            cusp::array1d	<Scalar     ,       cusp::host_memory>   	mu0(M,0);
            cusp::array1d	<Scalar     ,       cusp::host_memory>   	mu1(M,0);
            cusp::array1d	<Scalar     ,       cusp::host_memory>   	ranx(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	Phit(DD,0);   
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jn0(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jn1(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jm0(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jm1(DD,0);
			Scalar* pPhit=thrust::raw_pointer_cast(Phit.data());
			Scalar* pjn0 =thrust::raw_pointer_cast(jn0.data());
			Scalar* pjn1 =thrust::raw_pointer_cast(jn1.data());
			Scalar* pjm0 =thrust::raw_pointer_cast(jm0.data());
			Scalar* pjm1 =thrust::raw_pointer_cast(jm1.data());

            Scalar alpha=	 2.0f;
            Scalar beta =	-1.0f;
            Scalar zero =	 0.0f;
            Scalar one	=	 1.0f;
            Scalar mutemp;
			int		M0	=	(Floating)M/(Floating)2;
			for(Integer r=0;r<R;r++){
				//Define the Gaussian unnormalized random vector on host
				for(Integer i=0;i<DD;i++) 
				ranx[i]=0;	//ranx[i]=cusp::sqrt(-((Floating)2)*cusp::log((Scalar)rand()/(Scalar)RAND_MAX))*cusp::cos(((Floating)2*M_PI)*(Scalar)rand()/(Scalar)RAND_MAX);
							ranx[100]=1;
				//Pass this vector to the device name it jn0 and normalize it
				jn0=ranx;				
	            cublas::VectorTransform::normalize<Floating>(cublas_handle,DD,pjn0);//Normalize the vectors
				//Define the vector jm0=V|jn0>
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pjn0,&zero,pjm0);
				//Define the temporal vector vector |Phit>=V|jn0>
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pjn0,&zero,pPhit);
				//Calculate the first momento mu_00= <jn0|V^2|jn0> = <jm0|Phit>
                cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);//conjugate first vector
				mu[0*M+0]=mu[0*M+0]+ ((Floating)0.25)*mutemp;
				mu0[0]	 = mutemp;
				//Define the vector |jm1>= H|jm0>
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pjm0,&zero,pjm1);
				//Calculate the first momento mu_10= <jn0|VHV|jn0> = <jm1|Phit>
                cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhit,&mutemp);
				mu[1*M+0]=mu[1*M+0]+ ((Floating)0.5)*mutemp;
				mu0[1]	 = mutemp;
				for(int m=2;m<M;m++){
					//Define the vector |jm0>= 2H|jm1>-|jm0>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
					//Calculate the moments mu_{2,0},mu_{3,0},...,mu_{m,0}
					// as mu_{m,0}= <jn0|V T_m(H) V |jn0> = <jm0|Phit>
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
					mu[m*M+0]=mu[m*M+0]+((Floating)0.5)*mutemp;
					mu0[m]	 = mutemp;
					cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
					}
				//Define |jn1>=H|jn0> and |Phit>=V|jn1>
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH ,pjn0,&zero,pjn1);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pjn1,&zero,pPhit);
				//Define |jm0>=V|jn0> 
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pjn0,&zero,pjm0);
				//Calculate the moment mu_{0,1}=<jm0|Phit> 
				cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
				mu[0*M+1]	=mu[0*M+1]+ ((Floating)0.5)*mutemp;
				mu1[0]		= mutemp;
				//Define |jm1>=H|jm0>
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH ,pjm0,&zero,pjm1);
				//Calculate the moment mu_{1,1}=<jm0|Phit> 
				cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhit,&mutemp);				
				mu[1*M+1]	=mu[1*M+1]+ mutemp;
				mu1[1]		= mutemp;
				for(int m=2;m<M;m++){
					//Define the vector |jm0>= 2H|jm1>-|jm0>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
					//Calculate the moments mu_{2,1},mu_{3,1},...,mu_{m,1}
					// as mu_{m,1}= <jn0|V T_m(H) V H |jn0> = <jm0|Phit>
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
					mu[m*M+1]=mu[m*M+1]+mutemp;
					mu1[m]	 = mutemp;
					cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
					}
				//Define the vector |jm0>=V|jn1>
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GV,pjn1,&zero,pjm0);
				//Calculate the moment mu_{0,2} by using
				//mu_{0,2}=2 <Phi|T_1(H) V^2 T_1(H)|Phi>-\mu_{0,0}
				//mu_{0,2}=2 <jm0|Phit>-mu_{0,0}
				cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
				mu[0*M+2]=mu[0*M+2]+((Floating)0.5)*(((Floating)2)*mutemp-mu0[0]);
				//Define the vector |jm1>= H|jm0>
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pjm0,&zero,pjm1);
				//We calculate the moment mu_{1,2} by using
				//mu_{1,2}=2 <Phi|T_1(H)V HV_yT (H)|Phi>-\mu_{1,0}
				//mu_{1,2}=2 <jm1|Phit>-mu_{1,0}
				cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhit,&mutemp);				
				mu[1*M+2]=mu[1*M+2]+ (((Floating)2)*mutemp-mu0[1]);
				for(int m=2;m<M;m++){
					//Define the vector |jm0>= 2H|jm1>-|jm0>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
					//We calculate the moment mu_{m,2} by using
					//mu_{m,2}=2 <Phi|T_1(H)V T_m(H)V T_1(H)|Phi>-\mu_{m,0}
					//mu_{m,2}=2 <jm0|Phit>-mu_{m,0}
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);				
					mu[m*M+2]=mu[m*M+2]+(((Floating)2)*mutemp-mu0[m]);
					cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
					}
				for(int n=2;n<M0;n++){
					//Define the vector |jn0>= 2H|jn1>-|jn0>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH ,pjn1,&beta,pjn0);
					//Define the vector |Phit>= V|jn1>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GV,pjn1,&zero,pPhit);
					//Define the vector |jm0>= V|jn0>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GV,pjn0,&zero,pjm0);
					//We calculate the moment mu_{0,2n-1}
					//mu_{0,2n-1}=2 <Phi|T_n(H)V V T_{n-1}(H)|Phi>-\mu_{0,1}
					//mu_{0,2n-1}=2 <jm0|Phit>-mu_{0,1}
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
					mu[0*M+2*n-1]=mu[0*M+2*n-1]+((Floating)0.5)*(((Floating)2)*mutemp-mu1[0]);
					//Define the vector |jm1>= H|jm0>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GH ,pjm0,&zero,pjm1);
					//We calculate the moment mu_{1,2n-1}
					//mu_{1,2n-1}=2 <Phi|T_n(H)V  H V T_{n-1}(H)|Phi>-\mu_{1,1}
					//mu_{1,2n-1}=2 <jm1|Phit>-mu_{1,1}
					cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhit,&mutemp);
					mu[1*M+2*n-1]=mu[1*M+2*n-1]+(((Floating)2)*mutemp-mu1[1]);
					for(int m=2;m<M;m++){
						//Define the vector |jm0>= 2H|jm1>-|jm0>
						cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
						//We calculate the moment mu_{m,2n-1}
						//mu_{m,2n-1}=2 <Phi|T_n(H)V T_m(H) V T_{n-1}(H)|Phi>-\mu_{m,1}
						//mu_{m,2n-1}=2 <jm1|Phit>-mu_{m,1}
						cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
						mu[m*M+2*n-1]=mu[m*M+2*n-1]+(((Floating)2)*mutemp-mu1[m]);
						cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
					}
					//Define the vector |Phit>= V|jn0>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GV,pjn0,&zero,pPhit);
					//Define the vector |jm0>= V|jn0>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GV,pjn0,&zero,pjm0);
					//We calculate the moment mu_{0,2n}
					//mu_{0,2n}=2 <Phi|T_n(H)V^2 T_n(H)|Phi>-\mu_{0,0}
					//mu_{0,2n}=2 <jm0|Phit>-mu_{0,0}
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
					mu[0*M+2*n  ]=mu[0*M+2*n]+((Floating)0.5)*(((Floating)2)*mutemp-mu0[0]);
					//Define the vector |jm1>= H|jm0>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GH ,pjm0,&zero,pjm1);
					//We calculate the moment mu_{1,2n}
					//mu_{1,2n}=2 <Phi|T_n(H)V H V T_n(H)|Phi>-\mu_{1,0}
					//mu_{1,2n}=2 <jm1|Phit>-mu_{1,0}
					cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhit,&mutemp);
					mu[1*M+2*n]=mu[1*M+2*n]+(((Floating)2)*mutemp-mu0[1]);
					for(int m=2;m<M;m++){
						//Define the vector |jm0>= 2H|jm1>-|jm0>
						cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
						//We calculate the moment mu_{m,2n}
						//mu_{m,2n}=2 <Phi|T_2(H)V H V_yT_1(H)|Phi>-\mu_{m,0}
						//mu_{m,2n}=2 <jm1|Phit>-mu_{m,0}
						cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
						mu[m*M+2*n]=mu[m*M+2*n]+(((Floating)2)*mutemp-mu0[m]);
						cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
					}
				cublas::VectorTransform::swap(cublas_handle,DD,pjn0,pjn1);
				}
			}
        
        }

	template <typename Matrix, typename FloatType>
        void DOS( Matrix& H,typename Matrix::index_type M,typename Matrix::index_type R, const FloatType Emin,const  FloatType Emax,const FloatType alpha,const std::string outputname, int NE ){
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
			chebyshev::sum::dos(mu,Emin,Emax,alpha,(Floating)NE,outputname);   
        }

	template <typename Matrix,typename FloatType>
        void SIGMA( Matrix& H, Matrix& Vx,Matrix& Vy,int M, int R, const FloatType Emin,const  FloatType Emax,const FloatType alpha,const std::string outputnamexx,const std::string outputnamexy,int NE){ //This method will calculated the chebyshev moment using random vector method.
            typedef typename Matrix::value_type Scalar;
            typedef typename Matrix::index_type Integer;
            typedef typename Matrix::value_type::value_type Floating;

			
			std::string momentxx_out=outputnamexx;
			std::string momentxy_out=outputnamexy;

			size_t pos;
			pos = momentxx_out.find("data/");
			momentxx_out.replace(pos, std::string("data/").length(), "moments/");
			pos = momentxx_out.find(".dat");
			momentxx_out.replace(pos, std::string(".dat").length(), ".momXX");
			std::ofstream momentxx_file( momentxx_out.c_str() );

			pos = momentxy_out.find("data/");
			momentxy_out.replace(pos, std::string("data/").length(), "moments/");
			pos = momentxy_out.find(".dat");
			momentxy_out.replace(pos, std::string(".dat").length(), ".momXY");
			std::ofstream momentxy_file( momentxy_out.c_str() );

			std::cout<<"Storing the kpm moments in:"<<momentxx_out;

			cusp::array1d	<Scalar,cusp::host_memory> muxx_h(M*M,0);
			cusp::array1d	<Scalar,cusp::host_memory> muxy_h(M*M,0);
			cusp::array1d	<Scalar,cusp::device_memory> mu;

			std::cout<<"Calculating the moments"<<std::endl;
			chebyshev::random::GreenFun2D_moments(H,Vx,Vy,muxx_h,muxy_h,M,R);
			mu=muxx_h;

			for(int m=0;m< M; m++)
				for(int n=0;n< M; n++)
				{
					Floating scal=1./R; 
					
					momentxx_file<<m<<" "<<n<<" "<<scal*muxx_h[ m*M+n].real()<<" "<<scal*muxx_h[ m*M+n].imag()<<std::endl;
					momentxy_file<<m<<" "<<n<<" "<<scal*muxy_h[ m*M+n].real()<<" "<<scal*muxy_h[ m*M+n].imag()<<std::endl;
				}

			momentxx_file.close();
			momentxy_file.close();			
/*
			chebyshev::kernel::jackson2D_kernel(mu,1.0f/(Floating)R);            
			chebyshev::sum::Spectral_DCConductivity(mu,Emin,Emax,alpha,(Floating)NE,outputnamexx);
			mu=muxy_h;

			std::cout<<"SIGMAXY"<<std::endl;
			for(int m=0;m< M; m++)
				for(int n=0;n< M; n++)
					std::cout<<m<<" "<<n<<" "<<muxy_h[ m*M+n].real()<<" "<<muxy_h[ m*M+n].imag()<<std::endl;

			chebyshev::kernel::jackson2D_kernel(mu,1.0f/(Floating)R);            
			chebyshev::sum::Spectral_DCConductivity(mu,Emin,Emax,alpha,(Floating)NE,outputnamexy);
*/
			}




    }
    //Aqui se acaba chebyshev::random
namespace direct{
        
        template <typename Vector,typename IntType, typename FloatType>
        void PrintDOS(Vector& MDOS,Vector& TDOS, const IntType NE,const IntType R,const FloatType Emin, const FloatType Emax,const FloatType E0min, const FloatType E0max,const std::string outputname ){
            
            
            IntType k=0;
            FloatType dEinv=2/(Emax-Emin);
            FloatType dE   =(E0max-E0min)/(1.0f*NE);
            FloatType E,LDOS0, LDOS1 ;
            
            
            std::ofstream output_file;   output_file.open(outputname.c_str());
            
            for(E=E0min+2*dE;E<=E0max-2*dE;E=E+dE){
                LDOS0=          MDOS[k] *dEinv/((FloatType)R);
                LDOS1=cusp::exp(TDOS[k]/((FloatType)R))*dEinv;
                output_file<<E<<" "<<LDOS0<<" "<<LDOS1<<std::endl;
                k=k+1;
            }
            output_file.close();
            
            
        }
        
        
        template <typename Vector,typename FVector,typename IntType, typename FloatType>
        void ImGreen( Vector& mu, FVector& MDOS, FVector& TDOS,const IntType M, const IntType NE, const FloatType Emin,const  FloatType Emax,const FloatType E0min, const FloatType E0max){
            cusp::complex<FloatType> I(0.0f,1.0f);
            cusp::complex<FloatType> Gii;
            
            IntType m;
            IntType k=0;
            FloatType Enmin=(2*E0min-(Emax+Emin))/(Emax-Emin);
            FloatType Enmax=(2*E0max-(Emax+Emin))/(Emax-Emin);
            FloatType dEn  =(Enmax-Enmin)/(1.0f*NE);
            FloatType En, g;
            FloatType pi=3.14159265;
            FloatType alpha=2;
            FloatType LDOS;
            
            
            for(En=Enmin;En<Enmax;En=En+dEn){
                Gii=0.0; LDOS=0;
                for(m=0;m<M;m++){
                    g=(1/((FloatType)(M+1)))*(((FloatType)(M-m+1))*cusp::cos(pi*((FloatType)m)/((FloatType)(M+1)))+cusp::sin(pi*((FloatType)m)/((FloatType)(M+1)))*cusp::cos(pi/((FloatType)(M+1)))/cusp::sin(pi/((FloatType)(M+1))));
                    Gii=Gii+mu[m]*g*cusp::exp(I*((FloatType)m)*cusp::acos(En));
                }
                Gii=-alpha*I*Gii/pi;
                LDOS=-Gii.imag()/cusp::sqrt(1-En*En);
                MDOS[k]=MDOS[k]+LDOS;
                TDOS[k]=TDOS[k]+cusp::log(LDOS);
                //if(LDOS<0){std::cout<<"DOS NEGATIVA E="<<(2*En-(Emax+Emin))/(Emax-Emin)<<"DOS= "<<LDOS<<std::endl;}
                k=k+1;
            }
            
        }
        
        
        template <typename Matriz, typename FloatType>
        void DOS( Matriz& H,typename Matriz::index_type M,  typename Matriz::index_type Ns, const FloatType Emin,const FloatType Emax,const  FloatType E0min, const FloatType E0max, const std::string outputname, const std::string histname ){
            /********************************Algoritmo de Chebyshev y vectores aleatorios ***********************************/
            typedef typename Matriz::value_type Scalar;
            typedef typename Matriz::index_type Entero;
            Entero DD;
            FloatType pi=3.14159265;
            FloatType alpha=2;
            if(H.num_cols!=H.num_rows){ std::cout<<"ERROR: La matriz debe ser cuadrada, programa abortado"<<std::endl; exit(0);}	else {DD=H.num_cols;  }
            
            cusp::complex<FloatType> I(0.0,1.0);
            cusp::complex<FloatType> Gii;
            cusp::hyb_matrix<Entero, Scalar, cusp::host_memory>   H1=H;
            cusp::hyb_matrix<Entero, Scalar, cusp::device_memory> GH=H1;
            cusp::array1d	<Scalar     , 		cusp::device_memory> zero(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> j0(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> jn(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> jn1(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> jn2(DD,0);
            cusp::array1d	<Scalar     , 		cusp::host_memory> mu(M,0);
            
            const Entero NE=11213;
            FloatType MDOS[NE];
            FloatType TDOS[NE];
            for(int h=0;h<NE;h++) {MDOS[h]=0; TDOS[h]=0; }
            //FloatType    Enmin  =chebyshev::rescale(E0min,Emin,Emax);
            //FloatType    Enmax  =chebyshev::rescale(E0max,Emin,Emax);
            //FloatType    Endirac=chebyshev::rescale((FloatType)0,Emin,Emax);
            //FloatType   dEn     =(Enmax-Enmin)/((FloatType)1.0*NE);
            FloatType    Enmin  =0;
            FloatType    Enmax  =0;
            FloatType    Endirac=0;
            FloatType   dEn     =0;

            std::ofstream output_file;   output_file.open(histname.c_str());
            std::cout<<"Las energias para la distribucion son "<<std::endl;
            for(FloatType En=-300*dEn;En<=300*0.999*dEn;En=En+100*dEn)
                std::cout<<0.5*(0.999*(Emax-Emin)+(Emax+Emin))*En<<" ";
            std::cout<<std::endl;
            int l;
            for(int h=0;h<Ns;h++){
                cycletime(Ns);
                l=rand()%DD;
                jn=zero;
                jn1=zero;
                jn2=zero;
                jn[l]   =1.0f;
                mu[0]=0.5;
                multiply(GH,jn,jn1);
                mu[1]=jn1[l];
                for(int m=2;m<M;m++){
                    //printf("aqui estoy mundo\n");
                    multiply(GH,jn1,jn2);
                    cusp::VectorTransform::saxmy(alpha,jn2,jn);
                    mu[m]=jn2[l];
                    jn  = jn1;
                    jn1 = jn2;
                }
                ImGreen(mu,MDOS,TDOS,M,NE,Emin,Emax,E0min,E0max);
                for(FloatType En=-300*dEn;En<=300*0.999*dEn;En=En+100*dEn){
                    Gii=0.0;
                    for(int m=0;m<M;m++)
                        Gii=Gii+mu[m]*(1/((FloatType)(M+1)))*(((FloatType)(M-m+1))*cusp::cos(pi*((FloatType)m)/((FloatType)(M+1)))+cusp::sin(pi*((FloatType)m)/((FloatType)(M+1)))*cusp::cos(pi/((FloatType)(M+1)))/cusp::sin(pi/((FloatType)(M+1))))*(cusp::exp(I*((FloatType)m)*cusp::acos(En)));
                    Gii=-alpha*I*Gii/pi;
                    output_file<<-Gii.imag()/cusp::sqrt(1-En*En)<<" ";
                }
                output_file<<std::endl;
                cycletime(Ns);
            }
            output_file.close();
            PrintDOS(MDOS,TDOS,NE,Ns,Emin,Emax,E0min,E0max,outputname );
            /********************************Salida de los datos ***********************************/
        }
        
        
    }
    //Aqui se acaba chebyshev::direct
}
//Aqui se acaba chebyshev


template <typename Matriz,typename Dim, typename FloatType>
void LDOS( Matriz& H,Dim Nx, Dim Ny,typename Matriz::index_type M, const FloatType Emin,const  FloatType Emax, const FloatType Edirac,int W ,int L,const std::string outputname){
    /********************************Algoritmo de Chebyshev y vectores aleatorios ***********************************/
    typedef typename Matriz::value_type Scalar;
    typedef typename Matriz::index_type Entero;
	typedef typename Matriz::value_type::value_type Floating;

    Entero m,h,DD,D;
    const Entero NE=11213;
    FloatType E0min=-1.1;
    FloatType E0max= 1.1;
    FloatType Enmin,Enmax,Endirac,x,y;
    
    if((Emax-Emin)>0.00001){
        Enmin=(2*E0min-(Emax+Emin))/(Emax-Emin);
        Enmax=(2*E0max-(Emax+Emin))/(Emax-Emin);
        Endirac=(2*Edirac-(Emax+Emin))/(Emax-Emin);
    }
    else{
        Enmin  =(E0min)/Emax;
        Enmax  =(E0max)/Emax;
        Endirac=(Edirac)/(Emax-Emin);
    }
    
    FloatType dEn  =(Enmax-Enmin)/(1.0*NE);
    FloatType E00  =(Enmax+Enmin)/2;
    
    //FloatType lambda=30; FloatType delta=lambda/((float)(M));
    FloatType pi=3.14159265;
    FloatType alpha=2;
    FloatType En;
    
	DD=H.num_cols; D=DD/2; 
    cusp::complex<FloatType> I(0.0,1.0);
    cusp::complex<FloatType> Gii;
    cusp::hyb_matrix<Entero, Scalar, cusp::host_memory> H1=H;
    cusp::hyb_matrix<Entero, Scalar, cusp::device_memory> GH=H1;
    cusp::array1d	<Scalar     , 		cusp::device_memory> zero(DD,0);
    cusp::array1d	<Scalar     , 		cusp::device_memory> j0(DD,0);
    cusp::array1d	<Scalar     , 		cusp::device_memory> jn(DD,0);
    cusp::array1d	<Scalar     , 		cusp::device_memory> jn1(DD,0);
    cusp::array1d	<Scalar     , 		cusp::device_memory> jn2(DD,0);
    cusp::array1d	<Scalar     , 		cusp::host_memory> mu(M,0);
    
    std::ofstream output_file;   output_file.open(outputname.c_str());
  	for(int i=0;i<L;i++)for(int j=0;j<W;j++)for(int z=0;z<2;z++){
        h=((i+Nx)%Nx)*Ny+j+D*z;      
        x=(sqrt(3)*(Floating)(i+0.5*((Floating)(j-z))));
        //if(x>=((L-1)*sqrt(3))) x=x-1.5*L*sqrt(3);      
        y=0.5*(3*((Floating)j)-((Floating)z));
        cycletime(2*L*W);
        jn=zero;
        jn[h]   =1.0f;
        mu[0]=0.5;
        multiply(GH,jn,jn1);
        mu[1]=jn1[h];
        multiply(GH,jn1,jn2);
        cusp::VectorTransform::saxmy(alpha,jn2,jn);
        mu[2]=jn2[h];
        jn  = jn1;
        jn1 = jn2;
        for(m=3;m<M;m++){
            multiply(GH,jn1,jn2);
            cusp::VectorTransform::saxmy(alpha,jn2,jn);
            mu[m]=jn2[h];
            jn  = jn1;
            jn1 = jn2;
        }
        En=Endirac;
        Gii=0.0;
        for(m=0;m<M;m++)
            Gii=Gii+mu[m]*(1/((FloatType)(M+1)))*(((FloatType)(M-m+1))*cusp::cos(pi*((FloatType)m)/((FloatType)(M+1)))+cusp::sin(pi*((FloatType)m)/((FloatType)(M+1)))*cusp::cos(pi/((FloatType)(M+1)))/cusp::sin(pi/((FloatType)(M+1))))*(cusp::exp(I*((FloatType)m)*cusp::acos(En)));
        Gii=-alpha*I*Gii/pi;
        output_file<<x<<" "<<y<<" "<<-Gii.imag()/cusp::sqrt(1-En*En)<<std::endl;
        cycletime(2*W*L);
		}
    output_file.close();
    /********************************Salida de los datos ***********************************/

}

