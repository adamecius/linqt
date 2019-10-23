



#ifndef GRAPHENE_CUSP_H
#define GRAPHENE_CUSP_H
#ifndef UTILIDADES_H
#include "utilidades.h"
#define UTILIDADES_H
#endif
#endif




/****This library constain the definitions for some usual lattices*****/

namespace graphene{


  
    template <typename Dim,typename Matrix>	void lattice(Dim Nx, Dim Ny, Matrix& H)
	{
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.

                    //    V1=(i,j+1)_B
                    //                /
                    //  V0=(i,j)_A __/
                    //               \
                    //   V5=(i+1,j)_B \
                    
                    //here k is related to i position of the H_ij matrix, while
                    //h is related to the jm in this way h[0]+k0*v is the actual position
                    //on the sparse array. in this case we assume
                    //FIRST NEIGHBORHS
                    //         B                A
                    //____________________________________
                    //
                    //      \ h=2             h=3
                    //   h=0 \ __ h=1           /
                    //       /            h=1__/h=0
                    //      / h=3              \
                    //                      h=2 \

        //We determine the number of diagonals  v, which is related to the number of neighborns
        Integer D,h,k0,k1;      
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
        D=H.num_cols/(2);
		Scalar E0=(double) 0;
		Scalar t1=(double)-1;
		Scalar t2=(double)0;
		Scalar t3=(double)0;
		std::cout<<"Nx="<<Nx<<" Ny="<<Ny<<std::endl;
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix

	//ZERO NEIGHBOR
			/*Onsite (i,j)B-->(i,j)B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=0;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )B	<----> h=0
			H.column_indices[h+k0*v]= k0;								//Site: (i  ,j  )B	<----> h=0
			H.values[h+k0*v]		= 0;
			/*Onsite (i,j)A-->(i,j)A*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			h=0;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k0;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= 0;
	//FIRST NEIGHBOR
			/*Hoping (i,j)A-->(i,j)B*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=1;														
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t1;
			/*Hoping (i,j)B-->(i,j)A*/
			h=1;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )B				
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )A
			H.values[h+k1*v]		= t1;
			/*Hoping (i,j)A-->(i+1,j)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=2;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			H.column_indices[h+k0*v]= k1;								//Site: (i+1,j  )B
			H.values[h+k0*v]		= t1;
			/*Hoping (i+1,j)B-->(i,j)A*/
			h=2;
			H.row_indices[h+k1*v]	= k1;								//Site: (i+1,j  )B				
			H.column_indices[h+k1*v]= k0  ;								//Site: (i  ,j  )A
			H.values[h+k1*v]		= t1;
			/*Hoping (i,j)A-->(i,j+1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=3;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j+1)B
			H.values[h+k0*v]		= t1;
			/*Hoping (i,j+1)B-->(i,j)A*/
			h=3;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j+1)B				
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )A
			H.values[h+k1*v]		= t1;
		if(v==10){
	//SECOND NEIGHBOR
			/*Second (i,j)A<-->(i-1,j  )A*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			h=4;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)A<-->(i+1,j  )A*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
 			h=5;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)A<-->(i  ,j-1)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny  ;								
			h=6;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)A<-->(i  ,j+1)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny  ;								
			h=7;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)A<-->(i-1,j+1)A*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny  ;								
			h=8;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)A<-->(i+1,j-1)A*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny  ;								
			h=9;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)B<-->(i-1,j  )B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=4;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)B-->(i+1,j  )B*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=5;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)B-->(i  ,j-1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h=6;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)A-->(i  ,j+1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=7;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)A<-->(i-1,j+1)B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=8;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
			/*Second (i,j)A<-->(i+1,j-1)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h=9;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )A	<----> h=0
			H.values[h+k0*v]		= t2;
		}
		}
		
	}    
    template <typename Dim, typename Matrix>	void Magnetic_Field(Dim Nx, Dim Ny, Matrix& H,typename Matrix::index_type nPHI)
	{

		//In order to add a magnetic field to the tight-binding graphene hamiltonian, we need to define
		//the Peierls phase phi=+- pi (x/a) Phi/Phi0, 	and the magnetic flux MF= Phi/Phi0=Ba²sqrt(3) e/(2h)
			
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;						//This one will be use to counting and dimension parameters
        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.

		cusp::complex<Floating> I(0.0f,1.0f);      
        //We determine the number of diagonals  v, which is related to the number of neighborns
        Integer D,h,k0,k1;      
		//Test if spin is used or not;
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
        D=Nx*Ny;
		//double phi=(((double)2)*((double)nPHI))/((double)Ny);
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){


		Scalar PhiA[3]	={cusp::exp(I*(((Floating)2)*((Floating)M_PI))*((Floating)(j*nPHI)/(Floating)Ny)),cusp::exp(-I*(((Floating)2)*((Floating)M_PI))*((Floating)(j*nPHI)/(Floating)Ny)),1};

//		if(j==0&&i==23)
//			std::cout<<"Phi= "<<PhiA[0]<<" "<<PhiA[1]<<" "<<PhiA[2]<<std::endl;
//		if(j==Ny-1&&i==23)
//			std::cout<<"Phi= "<<PhiA[0]<<" "<<PhiA[1]<<" "<<PhiA[2]<<std::endl;
			
			//We start Filling the hoppings of the Matrix
			/*Hoping (i,j)A-->(i,j)B*/						
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;				
			k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=1;														
			H.values[h+k0*v]		= H.values[h+k0*v]*PhiA[0];
			/*Hoping (i,j)B-->(i,j)A*/
			h=1;
			H.values[h+k1*v]		= cusp::conj(H.values[h+k0*v]);
//				/*Hoping (i,j)A-->(i+1,j)B*/
//			if(i<Nx-1){
				k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				h=2;
				H.values[h+k0*v]		= H.values[h+k0*v]*PhiA[1];
				/*Hoping (i+1,j)B-->(i,j)A*/
				h=2;
				H.values[h+k1*v]		= cusp::conj(H.values[h+k0*v]);
//				}
//				if(j<Ny-1){
				/*Hoping (i,j)A-->(i,j+1)B*/
				k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
				h=3;
				H.values[h+k0*v]		= H.values[h+k0*v]*PhiA[2];
				/*Hoping (i,j+1)B-->(i,j)A*/
				h=3;
				H.values[h+k1*v]		= cusp::conj(H.values[h+k0*v]);
//				}
		}
			
	}
	
    template <typename Dim,typename Matrix>	void positionX(Dim Nx, Dim Ny,Matrix& X){
	
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
        Integer D,k0;      
        //then we define both the total number of atoms in one triangular lattices D, and the total number of atoms DD.
        D=Nx*Ny;
        Floating a=1.0*sqrt(3);
        Floating b=0.5*sqrt(3);

		/**********FIRST NEIGHBORNS******************/
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			/*Onsite (i,j)B-->(i,j)B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			X.row_indices[k0]	= k0;								//Site: (i  ,j  )B	<----> h=0
			X.column_indices[k0]= k0;								//Site: (i  ,j  )B	<----> h=0
			X.values[k0]		= a*i+b*j;
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			X.row_indices[k0]	= k0;								//Site: (i  ,j  )B	<----> h=0
			X.column_indices[k0]= k0;								//Site: (i  ,j  )B	<----> h=0
			X.values[k0]		= a*i+b*(j-1);
		}
	
  }

    template <typename Dim,typename Matrix>	void positionY(Dim Nx, Dim Ny,Matrix& Y){
	
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
        Integer D,k0;      
        //then we define both the total number of atoms in one triangular lattices D, and the total number of atoms DD.
        D=Nx*Ny;
        Floating a= 1.5;
        Floating b=-0.5;
		/**********FIRST NEIGHBORNS******************/
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			/*Onsite (i,j)B-->(i,j)B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			Y.row_indices[k0]	= k0;								//Site: (i  ,j  )B	<----> h=0
			Y.column_indices[k0]= k0;								//Site: (i  ,j  )B	<----> h=0
			Y.values[k0]		= a*j;
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			Y.row_indices[k0]	= k0;								//Site: (i  ,j  )B	<----> h=0
			Y.column_indices[k0]= k0;								//Site: (i  ,j  )B	<----> h=0
			Y.values[k0]		= a*j+b;

		}
	
  }

    template <typename Dim,typename Matrix>	void velocityx(Dim Nx, Dim Ny,Matrix& H,Matrix& V){
	
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);      
					//    V1=(i,j+1)_B
                    //                /
                    //  V0=(i,j)_A __/
                    //               \
                    //   V5=(i+1,j)_B \

				//(i,j)A-->(i  ,j  )B  ==>  (x,y)--> (x-1  ,    y      )
				//(i,j)A-->(i  ,j+1)B  ==>  (x,y)--> (x+1/2,y+sqrt(3)/2)
				//(i,j)A-->(i+1,j  )B  ==>  (x,y)--> (x+1/2,y-sqrt(3)/2)
				//(i,j)B-->(i  ,j  )A  ==>  (x,y)--> (x+1  ,    y      )
				//(i,j)B-->(i  ,j-1)A  ==>  (x,y)--> (x-1/2,y+sqrt(3)/2)
				//(i,j)B-->(i-1,j  )A  ==>  (x,y)--> (x-1/2,y-sqrt(3)/2)
/*               _____________ 
*				|\     |      |
*				| \ AA |  AB  |
*				|  \   |      |
*				|   \  |      |    
*               |    \ |      |
*               |_____\|______|
*				|      |\     |
*				|      | \ BB |
*				|  BA  |  \   |
*				|      |   \  |    
*               |      |    \ |
*               |_____ |_____\|
*/ 
        Integer D,h,k0,k1;      
   		//Test if spin is used or not;
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
        //Here we check if the matrix is squared. Ultimately this should be generalized
        if(H.num_cols!=H.num_rows){ std::cout<<"ERROR: The matrix should be an squared matrices."<<std::endl; exit(0);}
        //then we define both the total number of atoms in one triangular lattices D, and the total number of atoms DD.
        D=H.num_cols/(2);
   		Floating hbar=1.0;
   		Scalar* dxA = new Scalar[v];
		Scalar* dxB = new Scalar[v];
		/**********FIRST NEIGHBORNS******************/
		//dia[4]={(onsite)0,(NNN)0, 1, 0};
		//dja[4]={(onsite)0,(NNN)0, 0, 1};
		//dib[4]={(onsite)0,(NNN)0,-1, 0};
		//djb[4]={(onsite)0,(NNN)0, 0,-1};
		Scalar tempdxA[]={ (Floating)0      	,(Floating)-sqrt(3)/2,(Floating) sqrt(3)/2,(Floating) 0		  ,
							 -(Floating)sqrt(3)	, (Floating)sqrt(3)  ,(Floating)-sqrt(3)/2,(Floating) sqrt(3)/2,
							  -(Floating) sqrt(3)/2,(Floating) sqrt(3)/2};

		for(int i=0;i<v;i++)dxA[i]=tempdxA[i]/hbar; 		
		for(int i=0;i<H.num_entries;i++) V.values[i]=0;

		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			/*Onsite (i,j)B-->(i,j)B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=0;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )B	<----> h=0
			V.column_indices[h+k0*v]= k0;								//Site: (i  ,j  )B	<----> h=0
			V.values[h+k0*v]		= 0;
			/*Onsite (i,j)A-->(i,j)A*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			h=0;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			V.column_indices[h+k0*v]= k0;								//Site: (i  ,j  )A	<----> h=0
			V.values[h+k0*v]		= 0;
			/*Hoping (i,j)A-->(i,j)B*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=1;														
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j)B-->(i,j)A*/
			h=1;
			V.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )B				
			V.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )A
			V.values[h+k1*v]		= cusp::conj(V.values[h+k0*v]);
			/*Hoping (i,j)A-->(i+1,j)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=2;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j  )B
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i+1,j)B-->(i,j)A*/
			h=2;
			V.row_indices[h+k1*v]	= k1;								//Site: (i+1,j  )B				
			V.column_indices[h+k1*v]= k0  ;								//Site: (i  ,j  )A
			V.values[h+k1*v]		= cusp::conj(V.values[h+k0*v]);
			/*Hoping (i,j)A-->(i,j+1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=3;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j+1)B
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			h=3;
			V.row_indices[h+k1*v]	= k1;								//Site: (i  ,j+1)B				
			V.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )A
			V.values[h+k1*v]		= cusp::conj(V.values[h+k0*v]);
		if(v==10){
			//Second Neighborns A Lattice->Lattice
			/*Hoping (i,j)A-->(i-1,j)A*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			h=4;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i-1,j  )A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i+1,j)A*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			h=5;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j+1)A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j  )A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i,j-1)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
			h=6;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j-1)A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny;								
			h=7;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j+1)A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i,j+1)B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny;								
			h=8;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i-1,j+1)A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
			h=9;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j-1)A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i-1,j)B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=4;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i-1,j  )A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i+1,j)A*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=5;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j+1)A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j  )A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i,j-1)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h=6;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j-1)A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=7;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j+1)A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i,j+1)B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=8;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i-1,j+1)A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h=9;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j-1)A
			V.values[h+k0*v]		= I*dxA[h]*H.values[h+k0*v];
		}
		}
  }

    template <typename Dim, typename Matrix>	void velocityy(Dim Nx, Dim Ny,Matrix& H,Matrix& V){
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);     
					//    V1=(i,j+1)_B
                    //                /
                    //  V0=(i,j)_A __/
                    //               \
                    //   V5=(i+1,j)_B \

				//(i,j)A-->(i  ,j  )B  ==>  (x,y)--> (x-1  ,    y      )
				//(i,j)A-->(i  ,j+1)B  ==>  (x,y)--> (x+1/2,y+sqrt(3)/2)
				//(i,j)A-->(i+1,j  )B  ==>  (x,y)--> (x+1/2,y-sqrt(3)/2)
				//(i,j)B-->(i  ,j  )A  ==>  (x,y)--> (x+1  ,    y      )
				//(i,j)B-->(i  ,j-1)A  ==>  (x,y)--> (x-1/2,y+sqrt(3)/2)
				//(i,j)B-->(i-1,j  )A  ==>  (x,y)--> (x-1/2,y-sqrt(3)/2)
/*               _____________ 
*				|\     |      |
*				| \ AA |  AB  |
*				|  \   |      |
*				|   \  |      |    
*               |    \ |      |
*               |_____\|______|
*				|      |\     |
*				|      | \ BB |
*				|  BA  |  \   |
*				|      |   \  |    
*               |      |    \ |
*               |_____ |_____\|
*/ 
        Integer D,h,k0,k1;      
   		//Test if spin is used or not;
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
        //Here we check if the matrix is squared. Ultimately this should be generalized
        //then we define both the total number of atoms in one triangular lattices D, and the total number of atoms DD.
        D=H.num_cols/2;
   		Floating hbar=1.0;
   		Scalar* dyA = new Scalar[v];
		Scalar* dyB = new Scalar[v];
		/**********FIRST NEIGHBORNS******************/
		//dia[4]={(onsite)0,(NNN)0, 1, 0};
		//dja[4]={(onsite)0,(NNN)0, 0, 1};
		//dib[4]={(onsite)0,(NNN)0,-1, 0};
		//djb[4]={(onsite)0,(NNN)0, 0,-1};
		Scalar tempdyA[]={	(Floating)0,(Floating) -1/2	, (Floating)-1/2, 1			,
							(Floating)0,(Floating)0		,-(Floating) 3/2,(Floating)3/2	,(Floating)3/2,-(Floating)3/2};

		for(int i=0;i<v;i++)dyA[i]=tempdyA[i]/hbar; 
		for(int i=0;i<H.num_entries;i++) V.values[i]=0;

		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			/*Onsite (i,j)B-->(i,j)B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=0;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )B	<----> h=0
			V.column_indices[h+k0*v]= k0;								//Site: (i  ,j  )B	<----> h=0
			V.values[h+k0*v]		= 0;
			/*Onsite (i,j)A-->(i,j)A*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			h=0;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A	<----> h=0
			V.column_indices[h+k0*v]= k0;								//Site: (i  ,j  )A	<----> h=0
			V.values[h+k0*v]		= 0;
			/*Hoping (i,j)A-->(i,j)B*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=1;														
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j)B-->(i,j)A*/
			h=1;
			V.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )B				
			V.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )A
			V.values[h+k1*v]		= cusp::conj(V.values[h+k0*v]);
			/*Hoping (i,j)A-->(i+1,j)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=2;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j  )B
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i+1,j)B-->(i,j)A*/
			h=2;
			V.row_indices[h+k1*v]	= k1;								//Site: (i+1,j  )B				
			V.column_indices[h+k1*v]= k0  ;								//Site: (i  ,j  )A
			V.values[h+k1*v]		= cusp::conj(V.values[h+k0*v]);
			/*Hoping (i,j)A-->(i,j+1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=3;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j+1)B
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			h=3;
			V.row_indices[h+k1*v]	= k1;								//Site: (i  ,j+1)B				
			V.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )A
			V.values[h+k1*v]		= cusp::conj(V.values[h+k0*v]);
		if(v==10){
			//Second Neighborns A Lattice->Lattice
			/*Hoping (i,j)A-->(i-1,j)A*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			h=4;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i-1,j  )A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i+1,j)A*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			h=5;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j+1)A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j  )A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i,j-1)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
			h=6;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j-1)A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny;								
			h=7;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j+1)A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i,j+1)B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny;								
			h=8;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i-1,j+1)A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
			h=9;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j-1)A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i-1,j)B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=4;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i-1,j  )A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i+1,j)A*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=5;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j+1)A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j  )A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i,j-1)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h=6;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j-1)A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=7;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i  ,j+1)A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j)A-->(i,j+1)B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=8;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i-1,j+1)A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
			/*Hoping (i,j+1)B-->(i,j)A*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h=9;
			V.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			V.column_indices[h+k0*v]= k1;								//Site: (i+1,j-1)A
			V.values[h+k0*v]		= I*dyA[h]*H.values[h+k0*v];
		}
		}					
		
  }  
  
     template <typename Dim, typename Matrix, typename FloatType>
    void TSiteDis(Dim Nx,Dim Ny,Matrix& H, const FloatType p,const typename Matrix::value_type& U,const typename Matrix::value_type& tU)
    {
        
        typedef typename Matrix::value_type Scalar;
        typedef typename Matrix::index_type Integer;
        typedef typename Matrix::memory_space MemorySpace;
  		typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);   
        Integer D,v,k0,k1,h;
   		//Test if spin is used or not;
        //Here we check if the matrix is squared. Ultimately this should be generalized
        if(H.num_cols!=H.num_rows){ std::cout<<"ERROR: The matrix should be an squared matrices."<<std::endl; exit(0);}
        //then we define both the total number of atoms in one triangular lattices D, and the total number of atoms DD.
        else {D=H.num_cols/(2);}
   		v=(Integer)((double)H.num_entries/((double)H.num_cols));

        //dia[10]={0,0, 1, 0};
        //dja[10]={0,0, 0, 1};
        //dib[10]={0,0,-1, 0};
        //djb[10]={0,0, 0,-1};
        
        //Effective Hopping Matrix
        Scalar TA[4]={U,tU,tU,tU};
        Scalar TB[4]={U,tU,tU,tU};

		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			if((Floating)rand()/(Floating)RAND_MAX<p){
				/*Onsite (i,j)A-->(i,j)A*/
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
				h=0;
				H.values[h+k0*v]		=H.values[h+k0*v]+U;
				k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				h=1;
				H.values[h+k0*v]		=H.values[h+k0*v]+tU;
				H.values[h+k1*v]		=conj(H.values[h+k0*v]);

				k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				h=2;
				H.values[h+k0*v]		=H.values[h+k0*v]+tU;
				H.values[h+k1*v]		=conj(H.values[h+k0*v]);

				k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
				h=3;
				H.values[h+k0*v]		=H.values[h+k0*v]+tU;
				H.values[h+k1*v]		=conj(H.values[h+k0*v]);
				
				/*Hoping (i,j)A-->(i,j)B*/						
			}
			if((Floating)rand()/(Floating)RAND_MAX<p){
				/*Onsite (i,j)A-->(i,j)B*/
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				h=0;
				H.values[h+k0*v]		=H.values[h+k0*v]+U;
				
				k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
				h=1;
				H.values[h+k0*v]		=H.values[h+k0*v]+tU;
				H.values[h+k1*v]		=conj(H.values[h+k0*v]);

				k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
				h=2;
				H.values[h+k0*v]		=H.values[h+k0*v]+tU;
				H.values[h+k1*v]		=conj(H.values[h+k0*v]);

				k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
				h=3;
				H.values[h+k0*v]		=H.values[h+k0*v]+tU;
				H.values[h+k1*v]		=conj(H.values[h+k0*v]);

			}
		
		}

	}

     template <typename Dim, typename Matrix, typename FloatType>
    void LatticePotential(Dim Nx,Dim Ny,Matrix& H, const FloatType p,const typename Matrix::value_type& U)
    {
        
        typedef typename Matrix::value_type Scalar;
        typedef typename Matrix::index_type Integer;
        typedef typename Matrix::memory_space MemorySpace;
		typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);   
        Integer D,v,k0;
   		//Test if spin is used or not;
        //Here we check if the matrix is squared. Ultimately this should be generalized
        if(H.num_cols!=H.num_rows){ std::cout<<"ERROR: The matrix should be an squared matrices."<<std::endl; exit(0);}
        //then we define both the total number of atoms in one triangular lattices D, and the total number of atoms DD.
        else {D=H.num_cols/(2);}
   		v=(Integer)((double)H.num_entries/((double)H.num_cols));
		//dia[10]={0,0, 1, 0};
        //dja[10]={0,0, 0, 1};
        //dib[10]={0,0,-1, 0};

		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			if((Floating)rand()/(Floating)RAND_MAX<p){
				/*Onsite (i,j)A-->(i,j)A*/
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
				H.values[k0*v]		=H.values[k0*v]+((Floating)0.5)*U;
			}
			if((Floating)rand()/(Floating)RAND_MAX<p){
				/*Onsite (i,j)A-->(i,j)B*/
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				H.values[k0*v]		=H.values[k0*v]-((Floating)0.5)*U;
			}
		
		}

	}

     template <typename Dim, typename Matrix, typename FloatType>
    void SpinOrbitVacancies(Dim Nx,Dim Ny,Matrix& H, const FloatType p,FloatType lambda,FloatType VacReg)
    {
        
        typedef typename Matrix::value_type Scalar;
        typedef typename Matrix::index_type Integer;
        typedef typename Matrix::memory_space MemorySpace;
  		typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);   
        Integer D,v,k0,k1,h0,h1;
   		//Test if spin is used or not;
        //Here we check if the matrix is squared. Ultimately this should be generalized
        if(H.num_cols!=H.num_rows){ std::cout<<"ERROR: The matrix should be an squared matrices."<<std::endl; exit(0);}
        //then we define both the total number of atoms in one triangular lattices D, and the total number of atoms DD.
        else {D=H.num_cols/(2);}
   		v=(Integer)((double)H.num_entries/((double)H.num_cols));
        lambda=lambda/((Floating)6.0*sqrt(3));

        //dia[10]={0,0, 1, 0};
        //dja[10]={0,0, 0, 1};
        //dib[10]={0,0,-1, 0};
        //djb[10]={0,0, 0,-1};
        

		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			if((Floating)rand()/(Floating)RAND_MAX<p){
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
				/*Second (i,j)A<-->(i,j  )B*/
				k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				h0=1;
				H.values[h0+k0*v]		=(Floating)VacReg;
				H.values[h0+k1*v]		=(Floating)VacReg;
				/*Second (i,j)A<-->(i+1,j  )B*/
				k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				h0=2;
				H.values[h0+k0*v]		=(Floating)VacReg;
				H.values[h0+k1*v]		=(Floating)VacReg;
				/*Second (i,j)A<-->(i,j+1)B*/
				k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
				h0=3;
				H.values[h0+k0*v]		=(Floating)VacReg;
				H.values[h0+k1*v]		=(Floating)VacReg;
				/*Second (i,j)A<-->(i+1,j  )A*/
				k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
				h0=5;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
				h1=4;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
				/*Second (i,j)A<-->(i  ,j-1)A*/
				k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny  ;								
				h0=6;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
				h1=7;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
				/*Second (i,j)A<-->(i+1,j-1)A*/
				k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny  ;								
				h0=8;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
				h1=9;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);

			}
			if(rand()/(double)RAND_MAX<p){
				/*First (i,j)B*/
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				h0=0;
				/*First (i,j)B-->(i,j)A*/
				k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
				h0=1;
				H.values[h0+k0*v]		=(Floating)VacReg;
				H.values[h0+k1*v]		=(Floating)VacReg;
				/*First (i,j)B-->(i-1,j)A*/
				k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
				h0=2;
				H.values[h0+k0*v]		=(Floating)VacReg;
				H.values[h0+k1*v]		=(Floating)VacReg;
				/*First (i,j)B-->(i,j-1)A*/
				k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
				h0=3;
				H.values[h0+k0*v]		=(Floating)VacReg;
				H.values[h0+k1*v]		=(Floating)VacReg;
							/*SUBLATTICE B*/
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				/*Second (i,j)B-->(i+1,j  )B*/
				k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				h0=4;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
				h1=5;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
				/*Second (i,j)B-->(i  ,j-1)B*/
				k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
				h0=7;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
				h1=6;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
				/*Second (i,j)A<-->(i-1,j+1)B*/
				k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
				h0=9;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
				h1=8;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);

				
			}
		
		}

	}


     template <typename Dim, typename Matrix>
    void Anderson(Dim Nx,Dim Ny,Matrix& H, const typename Matrix::value_type& U)
    {
        
        typedef typename Matrix::value_type Scalar;
        typedef typename Matrix::index_type Integer;
        typedef typename Matrix::memory_space MemorySpace;
		typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);   
        Integer D=Nx*Ny;
        int k0;
		Floating WU;
        Integer v=(Integer)((Floating)H.num_entries/((Floating)H.num_cols));
		for(int i=0;i<Nx;i++)for(int j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			/*Onsite (i,j)A-->(i,j)A*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			WU=(Floating)rand()/((Floating)RAND_MAX);
			WU=(2.0*WU-1.0);
			H.values[k0*v]		=H.values[k0*v]+WU*U;
			/*Onsite (i,j)A-->(i,j)B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			WU=(Floating)rand()/((Floating)RAND_MAX);
			WU=(2.0*WU-1.0);
			H.values[k0*v]		=H.values[k0*v]+WU*U;		
		}

	}
    		
    template <typename Dim, typename Matrix, typename Floatype>	void HaldaneHoneyComb(Dim Nx, Dim Ny, Matrix& H,Floatype lambda)
	{

		//In order to add a magnetic field to the tight-binding graphene hamiltonian, we need to define
		//the Peierls phase phi=+- pi (x/a) Phi/Phi0, 	and the magnetic flux MF= Phi/Phi0=Ba²sqrt(3) e/(2h)
			
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
		typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);      
        //We determine the number of diagonals  v, which is related to the number of neighborns
        Integer D,h0,h1,k0,k1;      
		//Test if spin is used or not;
        Integer v=(Integer)((Floating)H.num_entries/((Floating)H.num_cols));
        D=Nx*Ny;
        lambda=lambda/((Floating)6.0*sqrt(3));
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
	//SECOND NEIGHBOR
			/*Second (i,j)A<-->(i-1,j  )A*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			/*Second (i,j)A<-->(i+1,j  )A*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
 			h0=5;
			H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
			h1=4;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*Second (i,j)A<-->(i  ,j-1)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny  ;								
			h0=6;
			H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
			h1=7;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*Second (i,j)A<-->(i+1,j-1)A*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny  ;								
			h0=8;
			H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
			h1=9;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*SUBLATTICE B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			/*Second (i,j)B-->(i+1,j  )B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h0=4;
			H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
			h1=5;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*Second (i,j)B-->(i  ,j-1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h0=7;
			H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
			h1=6;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*Second (i,j)A<-->(i-1,j+1)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h0=9;
			H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
			h1=8;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
		}
			
	}
	
    template <typename Dim, typename Matrix>	void PolarizedIntrinsecSO(Dim Nx, Dim Ny, Matrix& H,double lambda,double p)
	{

		//In order to add a magnetic field to the tight-binding graphene hamiltonian, we need to define
		//the Peierls phase phi=+- pi (x/a) Phi/Phi0, 	and the magnetic flux MF= Phi/Phi0=Ba²sqrt(3) e/(2h)
			
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
		typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);   
        //We determine the number of diagonals  v, which is related to the number of neighborns
        Integer D,h0,h1,k0,k1;      
		//Test if spin is used or not;
        Integer v=(Integer)((Floating)H.num_entries/((Floating)H.num_cols));
        D=Nx*Ny;
        lambda=lambda/10.3923;
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			
			if((Floating)rand()/(Floating)RAND_MAX<p){
		//SECOND NEIGHBOR
				/*Second (i,j)A<-->(i-1,j  )A*/
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
				/*Second (i,j)A<-->(i+1,j  )A*/
				k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
				h0=5;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
				h1=4;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
				/*Second (i,j)A<-->(i  ,j-1)A*/
				k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny  ;								
				h0=6;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
				h1=7;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
				/*Second (i,j)A<-->(i+1,j-1)A*/
				k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny  ;								
				h0=8;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+I*lambda;
				h1=9;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
				/*Second (i,j)B<-->(i-1,j  )B*/
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				/*Second (i,j)B-->(i+1,j  )B*/
				k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
				h0=5;
				H.values[h0+k0*v]		=H.values[h0+k0*v]-I*lambda;
				h1=4;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
				/*Second (i,j)B-->(i  ,j-1)B*/
				k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
				h0=6;
				H.values[h0+k0*v]		=H.values[h0+k0*v]-I*lambda;
				h1=7;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
				/*Second (i,j)A<-->(i+1,j-1)B*/
				k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
				h0=8;
				H.values[h0+k0*v]		=H.values[h0+k0*v]-I*lambda;
				h1=9;
				H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			}
		}
			
	}


    template <typename Dim, typename Matrix>	void KekuleHoneyComb(Dim Nx, Dim Ny, Matrix& H,double lambda)
	{

		//In order to add a magnetic field to the tight-binding graphene hamiltonian, we need to define
		//the Peierls phase phi=+- pi (x/a) Phi/Phi0, 	and the magnetic flux MF= Phi/Phi0=Ba²sqrt(3) e/(2h)
			
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
		typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);        
        //We determine the number of diagonals  v, which is related to the number of neighborns
        Integer D,h0,h1,k0,k1;      
		//Test if spin is used or not;
        Integer v=(Integer)((Floating)H.num_entries/((Floating)H.num_cols));
        D=Nx*Ny;
        lambda=lambda/10.3923;
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
	//SECOND NEIGHBOR
			/*Second (i,j)A<-->(i-1,j  )A*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
			/*Second (i,j)A<-->(i+1,j  )A*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny  ;								
 			h0=5;
			H.values[h0+k0*v]		=+I*lambda;
			h1=4;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*Second (i,j)A<-->(i  ,j-1)A*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny  ;								
			h0=6;
			H.values[h0+k0*v]		=+I*lambda;
			h1=7;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*Second (i,j)A<-->(i+1,j-1)A*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny  ;								
			h0=8;
			H.values[h0+k0*v]		=+I*lambda;
			h1=9;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*Second (i,j)B<-->(i-1,j  )B*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			/*Second (i,j)B-->(i+1,j  )B*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h0=5;
			H.values[h0+k0*v]		=-I*lambda;
			h1=4;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*Second (i,j)B-->(i  ,j-1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h0=6;
			H.values[h0+k0*v]		=-I*lambda;
			h1=7;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
			/*Second (i,j)A<-->(i+1,j-1)B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h0=8;
			H.values[h0+k0*v]		=-I*lambda;
			h1=9;
			H.values[h1+k1*v]		=conj(H.values[h0+k0*v]);
		}
			
	}

    //End os the space graphene::
}

template <typename Matrix>
void RefineSparse( Matrix& H)
{
    typedef typename Matrix::value_type Scalar;
    typedef typename Matrix::index_type Integer;
    typedef typename Matrix::memory_space MemorySpace;
    
	Integer k,kr,i,DD;
    Scalar zero =0.0;
    Integer v=(Integer)((double)H.num_entries/(double)H.num_cols);
	kr=0;
    if(H.num_cols!=H.num_rows){ std::cout<<"ERROR: La Matrix debe ser cuadrada, programa abortado"<<std::endl; exit(0);}	else {  DD=H.num_cols;}
	int nnz=0;for(i=0;i<DD*v;i++) if(H.values[i]!=zero) nnz=nnz+1;
    std::cout<<" La proyeccion de memoria a utilizar sin refinar es "<<sizeof(H.values[0])*DD*v*1.0/1000000000<<"GB";
    std::cout<<" luego de refinar es "<<sizeof(H.values[0])*nnz*1.0/1000000000<<"GB";
	std::cout<<" la enconomia de  memoria fue "<<sizeof(H.values[0])*(DD*v-nnz)*1.0/1000000000<<"GB"<<std::endl;
	Matrix Hr(DD,DD,nnz);
	for(k=0;k<DD*v;k++)
    if(H.values[k]!=zero){
        Hr.values[kr]        =H.values[k];
        Hr.column_indices[kr]=H.column_indices[k];
        Hr.row_indices[kr]   =H.row_indices[k];
        kr=kr+1;
    }
   	Hr.sort_by_row_and_column();
	H=Hr;
	
}

template <typename Matrix>
void RefineSparse( Matrix& H, int** CutMat)
{
    typedef typename Matrix::value_type Scalar;
    typedef typename Matrix::index_type Integer;
    typedef typename Matrix::memory_space MemorySpace;
	Integer k,kr,i,DD,v;
    Scalar zero =0.0;
    v=(Integer)((double)H.num_entries/(double)H.num_cols);
	kr=0;
    if(H.num_cols!=H.num_rows){ std::cout<<"ERROR: La Matrix debe ser cuadrada, programa abortado"<<std::endl; exit(0);}	else {  DD=H.num_cols;}
	int nnz=0;for(i=0;i<DD*v;i++) if(H.values[i]!=zero) nnz=nnz+1;
    std::cout<<" La proyeccion de memoria a utilizar sin refinar es "<<sizeof(H.values[0])*DD*v*1.0/1000000000<<"GB";
    std::cout<<" luego de refinar es "<<sizeof(H.values[0])*nnz*1.0/1000000000<<"GB";
	std::cout<<" la enconomia de  memoria fue "<<sizeof(H.values[0])*(DD*v-nnz)*1.0/1000000000<<"GB"<<std::endl;
	int rowtemp;
	//int columntemp;
	Matrix Hr(CutMat[0][0],CutMat[0][0],nnz);
    Scalar vacval=zero;
	rowtemp=H.row_indices[0];
    Hr.sort_by_row_and_column();
	for(k=0;k<v*DD;k=k+v){
        vacval=zero;
        for(int vnew=0;vnew<v;vnew++)
        vacval=vacval+H.values[k+vnew];
        if(vacval!=0) for(int vnew=0;vnew<v;vnew++)if(H.values[k+vnew]!=zero){
            Hr.values[kr]        =H.values[k+vnew];
            Hr.column_indices[kr]=H.column_indices[k+vnew];
            Hr.row_indices[kr]   =H.row_indices[k+vnew];
            kr=kr+1;
        }
    }
    for(int i=0;i<nnz;i++)
    for(int j=0;j<CutMat[0][0];j++){
        if(Hr.column_indices[i]==CutMat[0][j+1]) Hr.column_indices[i]=j;
        if(Hr.row_indices[i]==CutMat[0][j+1]) Hr.row_indices[i]=j;
    }
    
    //	Hr.sort_by_row();
	H=Hr;
	
}

/*******************************Finding Extrema****************************************/
/*******************************Finding Extrema****************************************/
template <typename Matrix,typename FloatType>	void FindingExtrema( Matrix& H,FloatType& Emin,  FloatType& Emax, int dT)
{
    
    typedef typename Matrix::value_type Scalar;
    typedef typename Matrix::index_type Integer;
    typedef typename Matrix::memory_space MemorySpace;
    cusp::complex<FloatType> ONE(1.0,0.0);
    std::cout<<"Calculando los extremos"<<std::endl;
    Integer m,i,it,DD,sample;
    if(H.num_cols!=H.num_rows){ std::cout<<"ERROR: La Matrix debe ser cuadrada, programa abortado"<<std::endl; exit(0);}        else { DD=H.num_cols;}
    
    FloatType dE,EM;
    FloatType Tolerance=0.0001;
    
    FloatType corrE=1.0+dT*Tolerance;
    cusp::array1d       <Scalar, cusp::host_memory>   Phi0(DD,0);
    cusp::array1d   <Scalar, cusp::device_memory> Psi0(DD,0);
    cusp::array1d   <Scalar, cusp::device_memory> Psif(DD,0);
    cusp::coo_matrix<Integer, Scalar, cusp::device_memory> GH;
    cusp::coo_matrix<Integer, Scalar, cusp::device_memory> Em(DD,DD,DD);
        
    GH = H;
    
    it=500; sample=100;
    double VarE, MedE, DesvAbsE, E[sample];
    int nit=0;
    int nmax=10000;
    VarE=10;
    for(i=0;i<DD;i++){ Phi0[i]=(FloatType)rand();}
    Psi0 =Phi0;
    cusp::VectorTransform::normalize(Psi0);
    
    //    std::cout<<"Calculando Emin"<<std::endl;
    while(cusp::abs(VarE)>Tolerance){
        nit=nit+1;
        for(m=1;m<it-sample;m++){
            multiply(GH,Psi0,Psif);
            Psi0=Psif;
            cusp::VectorTransform::normalize(Psi0);
        }
        
        for(m=it-sample;m<it;m++){
            multiply(GH,Psi0,Psif);
            E[m-it+sample]=(cusp::VectorTransform::dot(Psi0,Psif)*ONE).real();
            Psi0=Psif;
            cusp::VectorTransform::normalize(Psi0);
        }
        
        MedE=0; for(m=0;m<sample;m++) MedE=MedE+E[m];  MedE=MedE/sample;
        VarE=0; for(m=0;m<sample;m++) VarE=VarE+(E[m]-MedE)*(E[m]-MedE);  VarE=sqrt(VarE/sample);
        DesvAbsE=0; for(m=0;m<sample;m++) DesvAbsE=DesvAbsE+cusp::abs(E[m]-MedE);  DesvAbsE=DesvAbsE/sample;
        
        
        if(nit>nmax){ VarE=0; corrE=2*dT*Tolerance;}
        //      std::cout<<nit*it<<" "<<MedE<<" "<<VarE<<" "<<DesvAbsE<<" "<<cusp::abs(VarE)<<" "<<Tolerance<<std::endl;
    }
    
    Emin=MedE;
    
    nit=0;
    VarE=10;
    for(i=0;i<DD;i++){ Phi0[i]=(FloatType)rand();}
    Psi0 =Phi0;
    cusp::VectorTransform::normalize(Psi0);
    //std::cout<<"Calculando Emax"<<std::endl;
    while(cusp::abs(VarE)>Tolerance){
        nit=nit+1;
        for(m=1;m<it-sample;m++){
            multiply(GH,Psi0,Psif);
            cusp::VectorTransform::saxpy(-Emin,Psi0,Psif);
            cusp::VectorTransform::normalize(Psi0);
        }
        
        for(m=it-sample;m<it;m++){
            multiply(GH,Psi0,Psif);
            E[m-it+sample]=(cusp::VectorTransform::dot(Psi0,Psif)*ONE).real();
            cusp::VectorTransform::saxpy(-Emin,Psi0,Psif);
            cusp::VectorTransform::normalize(Psi0);
        }
        
        MedE=0; for(m=0;m<sample;m++) MedE=MedE+E[m];  MedE=MedE/sample;
        VarE=0; for(m=0;m<sample;m++) VarE=VarE+(E[m]-MedE)*(E[m]-MedE);  VarE=sqrt(VarE/sample);
        DesvAbsE=0; for(m=0;m<sample;m++) DesvAbsE=DesvAbsE+cusp::abs(E[m]-MedE);  DesvAbsE=DesvAbsE/sample;
        if(nit>nmax){ VarE=0; corrE=2*dT*Tolerance;}
        //  std::cout<<nit*it<<" "<<MedE<<" "<<VarE<<" "<<DesvAbsE<<" "<<cusp::abs(VarE)<<" "<<Tolerance<<std::endl;
    }
    
    Emax=MedE;
    FloatType temp; if(Emin>Emax){temp=Emax; Emax=Emin; Emin=temp;}
    std::cout<<"Energias del sistema Emin="<<Emin<<" Emax="<<Emax<<std::endl;
    dE=(Emax-Emin)*corrE*0.5;
    EM=(Emax+Emin)*0.5;
    Emax=dE+EM;
    Emin=EM-dE;
    std::cout<<"Energias del sistema Emin="<<Emin<<" Emax="<<Emax<<std::endl;
    for(i=0;i<DD;i++){
        Em.values[i]        =EM;
        Em.column_indices[i]=i;
        Em.row_indices[i]   =i;
    }
    cusp::subtract(GH, Em, GH);
    cusp::VectorTransform::matrixscale(1/dE, GH);
    H=GH;
    
}


template <typename Matrix,typename FloatType>	void FindingExtrema( Matrix& H,FloatType& Emin,  FloatType& Emax, FloatType Res)
{
    
    typedef typename Matrix::value_type Scalar;
    typedef typename Matrix::index_type Integer;
    typedef typename Matrix::memory_space MemorySpace;
 	cusp::complex<FloatType> I(0.0,1.0);
	std::cout<<"Calculando los extremos"<<std::endl;
    Integer DD;
    DD=H.num_cols;
    cusparseHandle_t	cusparse_handle;		//Handle of cusparse
    cusparseMatDescr_t	cusparse_descr;			//Cusparse's Matrix Descriptor  
    cusparseHybMat_t	GH;						//Cusparse's Hybrid Matrix
    cublasHandle_t		cublas_handle;			//Handle of cublas
    cusparseCreate			(&cusparse_handle);	//Initialization of the cusparse's handler
    cusparseCreateMatDescr	(&cusparse_descr);	//Initialization of the descriptor			
    cusparseCreateHybMat	(&GH);				//Initialization of the Hybrid Matrix
    cublasCreate			(&cublas_handle);	//Initialization of the cublas's handler  
    cusparse::convertX2Hyb(H,DD,cusparse_handle,cusparse_descr,GH);
    cusp::array1d	<Scalar     , 		cusp::device_memory> Psi0(DD,0);
    cusp::array1d	<Scalar     , 		cusp::device_memory> Psif(DD,0);
    cusp::array1d	<Scalar		,		cusp::host_memory>   Phi0(DD,0);
    Scalar* pPsi0=thrust::raw_pointer_cast(Psi0.data());
    Scalar* pPsif=thrust::raw_pointer_cast(Psif.data());    
	int nR=10;									//Number of random vectors we will use to calculate the mean bounds
	int nma=10000;									//Number of random vectors we will use to calculate the mean bounds
	int nit=0;
	FloatType Diff=1;							//Diference
	Scalar Et0,Et1;
	Scalar alpha =	-1.0f;
	Scalar zero  =	 0.0f;
	Scalar one	 =	 1.0f;


	for(int r=0;r<nR;r++){
		Phi0=Psi0;
		for(int i=0;i<DD;i++){ Phi0[i]= Phi0[i]+((FloatType)10)*Res*cusp::exp(I*((FloatType)rand()/((FloatType)RAND_MAX)));}
		Psi0 =Phi0;
		nit=0;
		while(Diff>Res&&Diff<2-Res){
			nit=nit+1;
			cublas::VectorTransform::normalize<FloatType>(cublas_handle,DD,pPsi0);//Here we normalize j0 and then we pass it to jn
			cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pPsi0,&zero,pPsif);
            cublas::VectorTransform::normalize<FloatType>(cublas_handle,DD,pPsif);//Here we normalize j0 and then we pass it to jn
			cublas::VectorTransform::swap(cublas_handle,DD,pPsif,pPsi0);
			cublas::VectorTransform::axpy(cublas_handle,DD,&alpha,pPsi0,pPsif);
			cublas::VectorTransform::nrm2(cublas_handle,DD,pPsif,&Diff);
			if(nit>=nma) break;
			}
		Diff=1;	
		}

	cusparse::Multiply(cusparse_handle,cusparse_descr,&one	,GH,pPsi0,&zero,pPsif);
	cublas::VectorTransform::dot(cublas_handle,DD,pPsi0,pPsif,&Et0);
	if(Et0.real()>0) Emax=Et0.real(); else if(Et0.real()<0) Emin=Et0.real();
	for(int i=0;i<DD;i++) Phi0[i]=0;
	Psi0=Phi0;
	Et0=-Et0;
	Diff=1;	
	for(int r=0;r<nR;r++){
		Phi0=Psi0;
		for(int i=0;i<DD;i++){ Phi0[i]= Phi0[i]+((FloatType)10)*Res*cusp::exp(I*((FloatType)rand()/((FloatType)RAND_MAX)));}
		Psi0 =Phi0;
		Psif =Phi0;
		nit=0;
		while(Diff>0.1*Res&&Diff<2-0.1*Res){
			nit=nit+1;
			cublas::VectorTransform::normalize<FloatType>(cublas_handle,DD,pPsi0);//Here we normalize j0 and then we pass it to jn
			cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pPsi0,&Et0,pPsif);
            cublas::VectorTransform::normalize<FloatType>(cublas_handle,DD,pPsif);//Here we normalize j0 and then we pass it to jn
			cublas::VectorTransform::swap(cublas_handle,DD,pPsif,pPsi0);
			cublas::VectorTransform::axpy(cublas_handle,DD,&alpha,pPsi0,pPsif);
			cublas::VectorTransform::nrm2(cublas_handle,DD,pPsif,&Diff);
			Psif =Psi0;
			if(nit>=nma) break;
			}
		Diff=1;	
		}	
		cusparse::Multiply(cusparse_handle,cusparse_descr,&one	,GH,pPsi0,&zero,pPsif);
		cublas::VectorTransform::dot(cublas_handle,DD,pPsi0,pPsif,&Et1);
		if(Et1.real()>0) Emax=Et1.real(); else if(Et1.real()<0) Emin=Et1.real();

		std::cout<<" Initial Energies Emax="<<Emax<<" Emin="<<Emin<<std::endl;


}

/*******************************Finding Extrema****************************************/
template <typename Matrix,typename FloatType>	void ExtremaPowerMethod( Matrix& H,FloatType& Emin,  FloatType& Emax, int dT)
{
    
    typedef typename Matrix::value_type Scalar;
    typedef typename Matrix::index_type Integer;
    typedef typename Matrix::memory_space MemorySpace;
    cusp::complex<FloatType> ONE(1.0,0.0);
    std::cout<<"Calculando los extremos"<<std::endl;
    Integer m,i,it,DD,sample;
    if(H.num_cols!=H.num_rows){ std::cout<<"ERROR: La Matrix debe ser cuadrada, programa abortado"<<std::endl; exit(0);}        else { DD=H.num_cols;}
    
    FloatType dE,EM;
    FloatType Tolerance=0.0001;
    dT       =400;
    
    FloatType corrE=1.0+dT*Tolerance;
    cusp::array1d       <Scalar, cusp::host_memory>   Phi0(DD,0);
    cusp::array1d   <Scalar, cusp::device_memory> Psi0(DD,0);
    cusp::array1d   <Scalar, cusp::device_memory> Psif(DD,0);
    cusp::coo_matrix<Integer, Scalar, cusp::device_memory> GH;
    cusp::coo_matrix<Integer, Scalar, cusp::device_memory> Em(DD,DD,DD);
    GH = H;
    
    it=500; sample=100;
    double VarE, MedE, DesvAbsE, E[sample];
    int nit=0;
    int nmax=10000;
    VarE=10;
    for(i=0;i<DD;i++){ Phi0[i]=(FloatType)rand();}
    Psi0 =Phi0;
    cusp::VectorTransform::normalize(Psi0);
    
    //    std::cout<<"Calculando Emin"<<std::endl;
    while(cusp::abs(VarE)>Tolerance){
        nit=nit+1;
        for(m=1;m<it-sample;m++){
            multiply(GH,Psi0,Psif);
            Psi0=Psif;
            cusp::VectorTransform::normalize(Psi0);
        }
        
        for(m=it-sample;m<it;m++){
            multiply(GH,Psi0,Psif);
            E[m-it+sample]=(cusp::VectorTransform::dot(Psi0,Psif)*ONE).real();
            Psi0=Psif;
            cusp::VectorTransform::normalize(Psi0);
        }
        
        MedE=0; for(m=0;m<sample;m++) MedE=MedE+E[m];  MedE=MedE/sample;
        VarE=0; for(m=0;m<sample;m++) VarE=VarE+(E[m]-MedE)*(E[m]-MedE);  VarE=sqrt(VarE/sample);
        DesvAbsE=0; for(m=0;m<sample;m++) DesvAbsE=DesvAbsE+cusp::abs(E[m]-MedE);  DesvAbsE=DesvAbsE/sample;
        
        
        if(nit>nmax){ VarE=0; corrE=2*dT*Tolerance;}
        //      std::cout<<nit*it<<" "<<MedE<<" "<<VarE<<" "<<DesvAbsE<<" "<<cusp::abs(VarE)<<" "<<Tolerance<<std::endl;
    }
    
    Emin=MedE;
    
    nit=0;
    VarE=10;
    for(i=0;i<DD;i++){ Phi0[i]=(FloatType)rand();}
    Psi0 =Phi0;
    cusp::VectorTransform::normalize(Psi0);
    //std::cout<<"Calculando Emax"<<std::endl;
    while(cusp::abs(VarE)>Tolerance){
        nit=nit+1;
        for(m=1;m<it-sample;m++){
            multiply(GH,Psi0,Psif);
            cusp::VectorTransform::saxpy(-Emin,Psi0,Psif);
            cusp::VectorTransform::normalize(Psi0);
        }
        
        for(m=it-sample;m<it;m++){
            multiply(GH,Psi0,Psif);
            E[m-it+sample]=(cusp::VectorTransform::dot(Psi0,Psif)*ONE).real();
            cusp::VectorTransform::saxpy(-Emin,Psi0,Psif);
            cusp::VectorTransform::normalize(Psi0);
        }
        
        MedE=0; for(m=0;m<sample;m++) MedE=MedE+E[m];  MedE=MedE/sample;
        VarE=0; for(m=0;m<sample;m++) VarE=VarE+(E[m]-MedE)*(E[m]-MedE);  VarE=sqrt(VarE/sample);
        DesvAbsE=0; for(m=0;m<sample;m++) DesvAbsE=DesvAbsE+cusp::abs(E[m]-MedE);  DesvAbsE=DesvAbsE/sample;
        if(nit>nmax){ VarE=0; corrE=2*dT*Tolerance;}
        //  std::cout<<nit*it<<" "<<MedE<<" "<<VarE<<" "<<DesvAbsE<<" "<<cusp::abs(VarE)<<" "<<Tolerance<<std::endl;
    } 
    Emax=MedE;
    FloatType temp; if(Emin>Emax){temp=Emax; Emax=Emin; Emin=temp;}
	std::cout<<"The boundaries of the Hamiltonian are Emin="<<Emin<<" "<<"Emax="<<Emax<<std::endl;
    
}


/*******************************Finding Extrema****************************************/
template <typename Matrix,typename FloatType>	void ExtremaExpMethod( Matrix& H,FloatType& Emin,  FloatType& Emax, FloatType& beta0)
{
    
            /********************************Algoritmo de Chebyshev y vectores aleatorios ***********************************/
            typedef typename Matrix::value_type Scalar;
            typedef typename Matrix::index_type Integer;
			typedef typename Scalar::value_type Floating;

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
            cusp::array1d	<Scalar     , 		cusp::device_memory	> Phi0(DD,0);
            cusp::array1d	<Scalar     , 		cusp::device_memory	> Phi1(DD,0);
  			cusp::array1d	<Scalar     ,       cusp::host_memory	> ranx(DD,0);
			Scalar* pPhi0=thrust::raw_pointer_cast(Phi0.data());
			Scalar* pPhi1=thrust::raw_pointer_cast(Phi1.data());
			Scalar temp;
            Scalar beta =	 beta0;
			int	   nmax	=	 3000;
            Scalar dbeta =	 beta/((Floating)nmax);
            Scalar zero =	 0.0f;
            Scalar one	=	 1.0f;
            //and then set the first moment which is always the same
			for(Integer i=0;i<DD;i++) ranx[i]=((Floating)rand())/((Floating)RAND_MAX);
			Phi0=ranx;
			cublas::VectorTransform::normalize(cublas_handle,DD,pPhi0);//Here we normalize j0 and then we pass it to jn
			Phi1=Phi0;									//We pass this vector to the first chebyshev vector jn
			for(int n=0;n<nmax;n++){
						cublas::VectorTransform::dot(cublas_handle,DD,pPhi0,pPhi1,&temp);
						cusparse::Multiply(cusparse_handle,cusparse_descr,&dbeta,GH,pPhi0,&one,pPhi1);
						Phi0=Phi1;
					}
			cublas::VectorTransform::normalize(cublas_handle,DD,pPhi0);//Here we normalize j0 and then we pass it to jn
			cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pPhi0,&zero,pPhi1);
			cublas::VectorTransform::dot(cublas_handle,DD,pPhi0,pPhi1,&temp);
			Emax=temp.real();

			dbeta =	-dbeta;
			Phi0=ranx;
			cublas::VectorTransform::normalize(cublas_handle,DD,pPhi0);//Here we normalize j0 and then we pass it to jn
			Phi1=Phi0;									//We pass this vector to the first chebyshev vector jn
			for(int n=0;n<nmax;n++){
					cublas::VectorTransform::dot(cublas_handle,DD,pPhi0,pPhi1,&temp);
					cusparse::Multiply(cusparse_handle,cusparse_descr,&dbeta,GH,pPhi0,&one,pPhi1);
					Phi0=Phi1;
				}
			cublas::VectorTransform::normalize(cublas_handle,DD,pPhi0);//Here we normalize j0 and then we pass it to jn
			cusparse::Multiply(cusparse_handle,cusparse_descr,&one,GH,pPhi0,&zero,pPhi1);
			cublas::VectorTransform::dot(cublas_handle,DD,pPhi0,pPhi1,&temp);
			Emin=temp.real();
			std::cout<<"The boundaries of the Hamiltonian are Emin="<<Emin<<" "<<"Emax="<<Emax<<std::endl;
			cusparseDestroy			(cusparse_handle);
			cusparseDestroyMatDescr	(cusparse_descr);
			cusparseDestroyHybMat	(GH);
			cublasDestroy			(cublas_handle);
    
}
