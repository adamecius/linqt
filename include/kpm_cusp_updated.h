#pragma once

#ifndef KPM_CUSP_H
#define KPM_CUSP_H
#endif

#include <thrust/functional.h>
#include "utilidades.h"
#include <thrust/inner_product.h>


namespace chebyshev{


/*******************************Finding Extrema****************************************/
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

namespace utilities{
	template <typename Matrix,typename GVector, typename Vector>		void delta_vector( Matrix& H, GVector& GX,Vector& g,int M){ //This method will calculated the chebyshev moment using random vector method.
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

		Scalar* pjn0=thrust::raw_pointer_cast(jn0.data());
		Scalar* pjn1=thrust::raw_pointer_cast(jn1.data());
		Scalar* pGX =thrust::raw_pointer_cast(GX.data());
		Scalar alpha=	 2.0f;
		Scalar beta =	-1.0f;
		Scalar zero =	 0.0f;
		Scalar one	=	 1.0f;

		//We pass the normalized vector to |jn0>
		jn0=GX;
		//Set  |GX>=0
		GX=jn1;
		//Set  |GX>=|GX>+g0|jn0>
		cublas::VectorTransform::axpy(cublas_handle,DD,&g[0],pjn0,pGX);					
		//Define |jn1>=1.0*H|jn0>+0.0*|jn1>
		cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pjn0,&zero,pjn1);
		//set   |GX>=|GX>+g1|jn1>= (g0+g1 H)|GX>
		cublas::VectorTransform::axpy(cublas_handle,DD,&g[1],pjn1,pGX);
		//Iterate it until M is reached
		for(Integer m=2;m<M;m++){
			// Redefine the vector |jn0>=2H|jn1>-|jn0> 
			cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjn1,&beta,pjn0);
			//Set |GX>=|GX>+gm|jn0>= (g0 + g1 H + gm Tm(H))|GX>
			cublas::VectorTransform::axpy(cublas_handle,DD,&g[m],pjn0,pGX);
			cublas::VectorTransform::swap(cublas_handle,DD,pjn0,pjn1);
			}
		//Pass back the vector to X
		cusparseDestroy			(cusparse_handle);
		cusparseDestroyMatDescr	(cusparse_descr);
		cusparseDestroyHybMat	(GH);
		cublasDestroy			(cublas_handle);
		}

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
		Floating pi=M_PI;

		int M=sqrt(mu.size());
		cusp::complex<Floating> zero=(Floating)0;	
		Scalar SIGMA;
		Spectral_DCConductivity_binary<Scalar>		delta_chev(a,M,a);							//Definimos la operacion binaria del kernel (mu,Index)
		thrust::plus<Scalar>	binary_op1;
		//Calculamos el kernel 
		std::ofstream output_file;   output_file.open(outputname.c_str());
		for(Floating En=-alpha;En<=alpha;En=En+dEn){
			a=2*pi/(pow((Emax-Emin)/(2*alpha),3)*pow(1-En*En,2));
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
				for(Integer i=0;i<DD;i++) HOSTphi[i]=cusp::sqrt(-((Floating)2)*cusp::log((Scalar)rand()/(Scalar)RAND_MAX))*cusp::cos((((Floating)2)*((Floating)M_PI))*(Scalar)rand()/(Scalar)RAND_MAX);			
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
				muxx[0*M+0]=muxx[0*M+0]+ ((Floating)0.25)*mutemp;
                cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);//conjugate first vector
				muxy[0*M+0]=muxy[0*M+0]+ ((Floating)0.25)*mutemp;
				cycletime((Floating)(M*M*R));
				// m=1     n=0//
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pjm0,&zero,pjm1);
				cycletime((Floating)(M*M*R));
                cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhitxx,&mutemp);
				muxx[1*M+0]=muxx[1*M+0]+ ((Floating)0.5)*mutemp;
                cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhitxy,&mutemp);
				muxy[1*M+0]=muxy[1*M+0]+ ((Floating)0.5)*mutemp;
				cycletime((Floating)(M*M*R));
				// all n's left  and m=0//
				for(int m=2;m<M;m++){
					cycletime((Floating)(M*M*R));
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxx,&mutemp);
					muxx[m*M+0]=muxx[m*M+0]+((Floating)0.5)*mutemp;
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);
					muxy[m*M+0]=muxy[m*M+0]+((Floating)0.5)*mutemp;
					cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
       				cycletime((Floating)(M*M*R));
					}
				//We now goes to the following m=1, and restart n=0;
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH ,pjn0,&zero,pjn1);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GVy,pjn1,&zero,pPhitxx);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GVx,pjn1,&zero,pPhitxy);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GVx,pPhi,&zero,pjm0);
				cycletime((Floating)(M*M*R));
				cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxx,&mutemp);
				muxx[0*M+1]=muxx[0*M+1]+ ((Floating)0.5)*mutemp;
				cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);
				muxy[0*M+1]=muxy[0*M+1]+ ((Floating)0.5)*mutemp;
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
					muxx[0*M+n]=muxx[0*M+n]+((Floating)0.5)*mutemp;
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhitxy,&mutemp);
					muxy[0*M+n]=muxy[0*M+n]+((Floating)0.5)*mutemp;
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

	for(int m=0;m<M;m++)
		for(int n=0;n<M;n++){
			muxy[m*M+n]=(muxy[m*M+n]+cusp::conj(muxy[n*M+m]))/(Floating)2;
			muxy[n*M+m]=cusp::conj(muxy[m*M+n]);
			muxx[m*M+n]=((muxx[m*M+n]+cusp::conj(muxx[n*M+m]))/(Floating)2).real();
			muxx[n*M+m]=muxx[m*M+n];
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
			cusparse::convertX2Hyb(V ,DD,cusparse_handle,cusparse_descr,GV);
            //Definitions of the device and hosts vectors
            cusp::array1d	<Scalar     ,       cusp::host_memory>   	HOSTphi(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	Phit(DD,0);   
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	Phi(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jn0(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jn1(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jm0(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory> 	jm1(DD,0);
			Scalar* pPhi =thrust::raw_pointer_cast(Phi.data());
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
			for(Integer r=0;r<R;r++){
				for(Integer i=0;i<DD;i++) HOSTphi[i]=cusp::sqrt(-((Floating)2)*cusp::log((Scalar)rand()/(Scalar)RAND_MAX))*cusp::cos((((Floating)2)*((Floating)M_PI))*(Scalar)rand()/(Scalar)RAND_MAX);			
				Phi=HOSTphi;
	            cublas::VectorTransform::normalize<Floating>(cublas_handle,DD,pPhi);//Normalize the vectors
				// m=0     n=0//
				//Ket
				jn0=Phi;
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pjn0,&zero,pPhit);
				//Bra
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pPhi,&zero,pjm0);
				cycletime((Floating)(M*M*R));
                cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);//conjugate first vector
				mu[0*M+0]=mu[0*M+0]+ ((Floating)0.25)*mutemp;
				cycletime((Floating)(M*M*R));
				// m=1     n=0//
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pjm0,&zero,pjm1);
				cycletime((Floating)(M*M*R));
                cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhit,&mutemp);
				mu[1*M+0]=mu[1*M+0]+ ((Floating)0.5)*mutemp;
				cycletime((Floating)(M*M*R));
				// all n's left  and m=0//
				for(int m=2;m<M;m++){
					cycletime((Floating)(M*M*R));
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
					mu[m*M+0]=mu[m*M+0]+((Floating)0.5)*mutemp;
					cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
       				cycletime((Floating)(M*M*R));
					}
				//We now goes to the following m=1, and restart n=0;
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH ,pjn0,&zero,pjn1);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV ,pjn1,&zero,pPhit);
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pPhi,&zero,pjm0);
				cycletime((Floating)(M*M*R));
				cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
				mu[0*M+1]=mu[0*M+1]+ ((Floating)0.5)*mutemp;
				cycletime((Floating)(M*M*R));
				// n=1     m=1//
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH ,pjm0,&zero,pjm1);
				cycletime((Floating)(M*M*R));
				cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhit,&mutemp);				
				mu[1*M+1]=mu[1*M+1]+ mutemp;
				cycletime((Floating)(M*M*R));
				// all m's left  and n=1//
				for(int m=2;m<M;m++){
					cycletime((Floating)(M*M*R));
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
					mu[m*M+1]=mu[m*M+1]+mutemp;
					cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
					cycletime((Floating)(M*M*R));
					}
				//We now iterate over all n's and restart m;
				for(int n=2;n<M;n++){
					cycletime((Floating)(M*M*R));
					cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH ,pjn1,&beta,pjn0);
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GV,pjn0,&zero,pPhit);
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GV,pPhi,&zero,pjm0);
					cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
					mu[0*M+n]=mu[0*M+n]+((Floating)0.5)*mutemp;
					cycletime((Floating)(M*M*R));
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one  ,GH ,pjm0,&zero,pjm1);
					cycletime((Floating)(M*M*R));
					cublas::VectorTransform::dot(cublas_handle,DD,pjm1,pPhit,&mutemp);
					mu[1*M+n]=mu[1*M+n]+ mutemp;
					cycletime((Floating)(M*M*R));
					for(int m=2;m<M;m++){
						cycletime((Floating)(M*M*R));
						cusparse::Multiply(cusparse_handle,cusparse_descr,&alpha,GH,pjm1,&beta,pjm0);
						cublas::VectorTransform::dot(cublas_handle,DD,pjm0,pPhit,&mutemp);
						mu[m*M+n]=mu[m*M+n]+mutemp;
						cycletime((Floating)(M*M*R));
						cublas::VectorTransform::swap(cublas_handle,DD,pjm0,pjm1);
					}					
				cublas::VectorTransform::swap(cublas_handle,DD,pjn0,pjn1);
				}
			}

	for(int m=0;m<M;m++)
		for(int n=0;n<M;n++){
			mu[m*M+n]=((mu[m*M+n]+cusp::conj(mu[n*M+m]))/(Floating)2).real();
			mu[n*M+m]=mu[m*M+n];
			}
        };

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

			cusp::array1d	<Scalar,cusp::host_memory> muxx_h(M*M,0);
			cusp::array1d	<Scalar,cusp::host_memory> muxy_h(M*M,0);
			cusp::array1d	<Scalar,cusp::device_memory> mu;

			chebyshev::random::GreenFun2D_moments(H,Vx,Vy,muxx_h,muxy_h,M,R);
			mu=muxx_h;
			chebyshev::kernel::jackson2D_kernel(mu,1.0f/(Floating)R);            
			chebyshev::sum::Spectral_DCConductivity(mu,Emin,Emax,alpha,(Floating)NE,outputnamexx);
			mu=muxy_h;
			chebyshev::kernel::jackson2D_kernel(mu,1.0f/(Floating)R);            
			chebyshev::sum::Spectral_DCConductivity(mu,Emin,Emax,alpha,(Floating)NE,outputnamexy);
			}

	template <typename Matrix,typename FloatType>
        void SIGMA( Matrix& H, Matrix& Vx,int M, int R, const FloatType Emin,const  FloatType Emax,const FloatType alpha,const std::string outputnamexx,int NE){ //This method will calculated the chebyshev moment using random vector method.
            typedef typename Matrix::value_type Scalar;
            typedef typename Matrix::index_type Integer;
            typedef typename Matrix::value_type::value_type Floating;

			cusp::array1d	<Scalar,cusp::host_memory> muxx_h(M*M,0);
			cusp::array1d	<Scalar,cusp::host_memory> muxy_h(M*M,0);
			cusp::array1d	<Scalar,cusp::device_memory> mu;

			chebyshev::random::GreenFun2D_moments(H,Vx,muxx_h,M,R);
			mu=muxx_h;
			chebyshev::kernel::jackson2D_kernel(mu,1.0f/(Floating)R);            
			chebyshev::sum::Spectral_DCConductivity(mu,Emin,Emax,alpha,(Floating)NE,outputnamexx);
			}


	template <typename Matrix,typename FloatType>
        void DosConductivity(Matrix& H ,Matrix& V,int M ,int R ,const FloatType Emin ,const FloatType Emax,const FloatType alpha, const std::string outputdos ,const std::string outputcond ,const  FloatType E0min ,const FloatType E0max,int NE){ //This method will calculated the chebyshev moment using random vector method.
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
			Floating hbar=	 0.65821;
			Scalar zero =	 0.0f;
			Scalar one	=	 1.0f;
 			//Definitions of integers variables
			Integer DD=H.num_cols;  			
			cusparse::convertX2Hyb(V,DD,cusparse_handle,cusparse_descr,GV);
            //Definitions of the device and hosts vectors
            cusp::array1d	<Scalar     , 		cusp::device_memory	> 	GP(DD,0);            
            cusp::array1d	<Scalar     , 		cusp::device_memory	> 	GX(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory	> 	GY(DD,0);
            cusp::array1d	<Scalar     , 		cusp::host_memory	> 	Psi(DD,0);
            cusp::array1d	<Scalar     , 		cusp::host_memory	> 	g(M ,0);
            cusp::array1d	<Floating   , 		cusp::host_memory	> 	Ener(NE ,0);
            cusp::array1d	<Scalar     , 		cusp::host_memory	> 	Dos (NE ,0);
            cusp::array1d	<Scalar     , 		cusp::host_memory	> 	Cond(NE ,0);

			Scalar* pGP =thrust::raw_pointer_cast(GP.data());
			Scalar* pGX =thrust::raw_pointer_cast(GX.data());
			Scalar* pGY =thrust::raw_pointer_cast(GY.data());
			std::ofstream output_dos ;   output_dos.open( outputdos.c_str());
			std::ofstream output_cond;   output_cond.open(outputcond.c_str());
			//Some necessary variables
			Floating dE=(E0max-E0min)/(Floating)NE;
			Floating  E=E0min;
			Floating pi		=(Floating)M_PI;
			Floating gJackson;
			Floating RenormE=2*alpha/(Emax-Emin);
			Scalar tempcond	=0;
			Scalar tempdos	=0;

			//Sum over all the  Random vectors |r>, r=0,1,2,...
			for(Integer r=0;r<R;r++){
				cycletime(NE*R);
				//Define the Gaussian unnormalized random vector on host
				for(Integer i=0;i<DD;i++) Psi[i]=cusp::sqrt(-((Floating)2)*cusp::log((Floating)rand()/(Floating)RAND_MAX))*cusp::cos(((Floating)2*M_PI)*(Floating)rand()/(Floating)RAND_MAX);					
				//SUM OVER ALL THE ENERGY PARTITION VARIABLE E=Emin+n*dE
				E=E0min;
				for( int n=0; n<=NE;n++){
					Ener[n]=E;				
					Floating En		=alpha*(2*E-(Emax+Emin))/(Emax-Emin);
					E=E+dE;
					//Define a vector of the kernel and energy-constant factors for  for all moments
					// V=( g0 T_0(En)/pi*sqrt(1-En^2), 
					for(int m=0;m<M;m++){
						gJackson=((M-m+1)*cos(pi*m/(M+1))+sin(pi*m/(M+1))*cos(pi/(M+1))/sin(pi/(M+1)))/(M+1);
						g[m]=gJackson*RenormE*cusp::cos(((Floating)m)*acos(En))*2.0/(pi*cusp::sqrt(1-En*En));
						}
					g[0]=g[0]*0.5;
					//Define the vector |GX>=|Psi>
					GP=Psi;
					//Define the vector |GX>=|Psi>/ <Psi|Psi>
					cublas::VectorTransform::normalize<Floating>(cublas_handle,DD,pGP);//Normalize the vectors
					GX=GP;
					//Define the vector |GX>=delta(E-H)|Psi>
					chebyshev::utilities::delta_vector(H,GX,g,M);
					//Define the vector dos=<Psi|delta(E-H)|Psi>
					cublas::VectorTransform::dot(cublas_handle,DD,pGP,pGX,&tempdos);//conjugate first vector
					Dos[n]=Dos[n]+tempdos;
					//Define the vector |GY> delta(E-H)V|GY>
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pGP,&zero,pGY);
					chebyshev::utilities::delta_vector(H,GY,g,M);
					//Define the vector |GX>=V delta(E-H) |GX>
					GP=GX;
					cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pGP,&zero,pGX);
					//Calculate the Conductivity as   sigmaxx= <Y| V delta(E-H) V delta(E-H) |X>
					cublas::VectorTransform::dot(cublas_handle,DD,pGY,pGX,&tempcond);//conjugate first vector
					Cond[n]=Cond[n]+tempcond;
					cycletime(NE*R);
					}
				}

			//Print the DOS and Cond ,...
			for(int n=0;n<NE;n++){
				Dos[n] =Dos[n] /(((Floating)R));
				Cond[n]=Cond[n]*(hbar*hbar*pi)/((Floating)R);
				output_dos <<Ener[n]<<" "<< (Dos[n] ).real()<<std::endl; 
				output_cond<<Ener[n]<<" "<< (Cond[n]).real()<<std::endl;   
				}
			
			
			cusparseDestroy			(cusparse_handle);
			cusparseDestroyMatDescr	(cusparse_descr);
			cusparseDestroyHybMat	(GV);
			cublasDestroy			(cublas_handle);
			}

    }

namespace realspace{
		template <typename Matrix,typename FloatType>
        void LinearDosConductivity(int Nx, int Ny,Matrix& H ,Matrix& V,int M ,const FloatType Emin ,const FloatType Emax,const FloatType alpha,const std::string outputdos ,const std::string outputcond ,const  FloatType E0,FloatType Lx,FloatType Ly){ //This method will calculated the chebyshev moment using random vector method.
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
			Floating hbar=	 0.65821;
			Scalar zero =	 0.0f;
			Scalar one	=	 1.0f;
 			//Definitions of integers variables
			Integer D=Nx*Ny;
			Integer DD=2*D;						
			cusparse::convertX2Hyb(V,DD,cusparse_handle,cusparse_descr,GV);
            //Definitions of the device and hosts vectors
            cusp::array1d	<Scalar     , 		cusp::device_memory	> 	GP(DD,0);            
            cusp::array1d	<Scalar     , 		cusp::device_memory	> 	GX(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory	> 	GY(DD,0);
            cusp::array1d	<Scalar     , 		cusp::host_memory	> 	Psi(DD,0);
            cusp::array1d	<Scalar     , 		cusp::host_memory	> 	Vzero(DD,0);
            cusp::array1d	<Scalar     , 		cusp::host_memory	> 	g(M ,0);
			Scalar* pGP =thrust::raw_pointer_cast(GP.data());
			Scalar* pGX =thrust::raw_pointer_cast(GX.data());
			Scalar* pGY =thrust::raw_pointer_cast(GY.data());
			std::ofstream output_dos ;   output_dos.open( outputdos.c_str());
			std::ofstream output_cond;   output_cond.open(outputcond.c_str());

			Floating  gJackson;
			Scalar cond,dos;
			Floating pi		=(Floating)M_PI;
			Floating RenormE=2*alpha/(Emax-Emin);
			Floating En		=alpha*(2*E0-(Emax+Emin))/(Emax-Emin);
			Integer h;
			for(int m=0;m<M;m++){
				gJackson=((M-m+1)*cos(pi*m/(M+1))+sin(pi*m/(M+1))*cos(pi/(M+1))/sin(pi/(M+1)))/(M+1);
				g[m]=gJackson*RenormE*cusp::cos(((Floating)m)*acos(En))*2.0/(pi*cusp::sqrt(1-En*En));
				}
			g[0]=g[0]*0.5;

			/*Because we want to calculate the DOS and Conductivity in a
			 * real space window, it is a good thing to choose the primitive
			 * vectors and use them to evaluate the actual position of the 
			 * atoms
			 */
			 Floating a1[]={ 1.0*sqrt(3), 0.0};
			 Floating a2[]={ 0.5*sqrt(3), 1.5};
			 Floating d1[]={-0.5*sqrt(3),-0.5};
			/***Center of the Disk initially fixed at the center of graphene. CAN BE EXTENDED***/
			int i0=(Integer)((Floating)Nx/2.0); /*Centered at the X direction*/
			int j0=(Integer)((Floating)Ny/2.0); /*Centered at the Y direction*/
			Floating rc[]={a1[0]*i0+a2[0]*j0,a1[1]*i0+a2[1]*j0};
			Floating r []={a1[0]*i0+a2[0]*j0,a1[1]*i0+a2[1]*j0};
		
			for(int i=0;i<Nx;i++)for(int j=0;j<Ny;j++)for(int z=0;z<2;z++){
			//*Calculate the r point associated with (i,j)_z indexes
			r[0]=a1[0]*i+a2[0]*j+d1[0]*z;
			r[1]=a1[1]*i+a2[1]*j+d1[1]*z;
			h=i*Ny+j+D*z;
			/*Check if the point belongs to the selected region
			 * In this case the intersecption between the two circles C1=(r0,Rmax) and C2=(r0,Rmin)
			 */
			cycletime(2*Nx*Ny);
			if((r[0]>rc[0]-Lx/2.0&&r[0]<rc[0]+Lx/2.0)&&(r[1]>rc[1]-Ly/2.0&&r[1]<rc[1]+Ly/2.0)){       
				cond		=0;	
				dos			=0; 
				//Define the Gaussian unnormalized random vector on host
				//In this case  the vector Psi= delta_h,0
				Psi		=Vzero;
				Psi[h]	=1.0;			
				GP=Psi;
				//Define the vector |GX>=|Psi>
				GX=GP;
				//Define the vector |GX>=delta(E-H)|Psi>
				chebyshev::utilities::delta_vector(H,GX,g,M);
				//Define the vector dos=<Psi|delta(E-H)|Psi>
                cublas::VectorTransform::dot(cublas_handle,DD,pGP,pGX,&dos);//conjugate first vector
				//Define the vector |GY> delta(E-H)V|GY>
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pGP,&zero,pGY);
				chebyshev::utilities::delta_vector(H,GY,g,M);
				//Define the vector |GX>=V delta(E-H) |GX>
				GP=GX;
				cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GV,pGP,&zero,pGX);
				//Calculate the Conductivity as   sigmaxx= <Y| V delta(E-H) V delta(E-H) |X>
                cublas::VectorTransform::dot(cublas_handle,DD,pGY,pGX,&cond);//conjugate first vector

				cond=cond*5.0*hbar*hbar*pi;
				output_dos <<" "<<r[0]<<" "<<r[1]<<" "<<dos .real()<<std::endl; 
				output_cond<<" "<<r[0]<<" "<<r[1]<<" "<<cond.real()<<std::endl;   
				}
				cycletime(2*Nx*Ny); 
			}
				
			cusparseDestroy			(cusparse_handle);
			cusparseDestroyMatDescr	(cusparse_descr);
			cusparseDestroyHybMat	(GV);
			cublasDestroy			(cublas_handle);
			}
	}

}
//Aqui se acaba chebyshev

