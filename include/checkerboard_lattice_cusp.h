



#ifndef GRAPHENE_CUSP_H
#define GRAPHENE_CUSP_H
#ifndef UTILIDADES_H
#include "utilidades.h"
#define UTILIDADES_H
#endif
#endif




/****This library constain the definitions for some usual lattices*****/

namespace square{

template <typename Dim,typename Matrix>	void lattice(Dim Nx, Dim Ny, Matrix& H)
	{
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.

                    //FIRST NEIGHBORHS

//	The CheckerBoard Lattice has two inequivalent sites * and º
//  and the lattice can be written in the following way
//               
//                       *       ºb2
//                       |       |
//                       |       |
//				 *_______º_______*______ºb3
//                       |b0     |a0
//                       |       |
//                       *b      ºb1
//                       |       |
//                       |       |
//				 *_______º_______*______º b
//                       |b	     |a
//                       |       |
//                       *       º b
	 
        //We determine the number of diagonals  v, which is related to the number of neighborns
		cusp::complex<Floating> I(0.0f,1.0f);      
        Integer h,k0,k1,D;
		double Phi=0.7853;
        D		=Nx*Ny;      
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
		Scalar E0=(double) 0;
		Scalar t1=(double)-1;
		Scalar t2=(double)-1/(sqrt(2)+2);
		Scalar t3=(double)-1/(2*sqrt(2)+2);

		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
	
	//ZERO NEIGHBOR  
			/*Onsite (i,j)b-->(i,j)b*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			h=0;	
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )B	<----> h=0
			H.column_indices[h+k0*v]= k0;								//Site: (i  ,j  )B	<----> h=0
			H.values[h+k0*v]		=-E0;
			/*Onsite (i,j)a-->(i,j)a*/
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			h=0;	
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )B	<----> h=0
			H.column_indices[h+k0*v]= k0;								//Site: (i  ,j  )B	<----> h=0
			H.values[h+k0*v]		= E0;
	//FIRST NEIGHBOR
			/*Hoping (i,j)a<-->(i+0,j+0)b*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			//Direct
			h=1;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t1*cusp::exp( I*Phi);
			//Inverse
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t1*cusp::exp(-I*Phi);
			/*Hoping (i,j)a<-->(i+1,j+0)b*/						
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			//Direct
			h=2;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t1*cusp::exp(-I*Phi);
			//Inverse
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t1*cusp::exp( I*Phi);
			/*Hoping (i,j)a<-->(i+0,j+1)b*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			//Direct
			h=3;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t1*cusp::exp(-I*Phi);
			//Inverse
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t1*cusp::exp( I*Phi);
			/*Hoping (i,j)a<-->(i+1,j+1)b*/						
			k1=((i+1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			//Direct
			h=4;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t1*cusp::exp( I*Phi);
			//Inverse
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t1*cusp::exp(-I*Phi);
	//SECOND NEIGHBOR
			/*Hoping (i,j)a<-->(i-1,j+0)*/						
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			//Direct
			h=5;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t2;
			//Inverse
			h=6;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t2;
			/*Hoping (i,j)a<-->(i+1,j+0)a*/						
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			//Direct
			h=6;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t2;
			//Inverse
			h=5;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t2;
			/*Hoping (i,j)a<-->(i+0,j-1)a*/						
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
			//Direct
			h=7;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t2;
			//Inverse
			h=8;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t2;
			/*Hoping (i,j)a<-->(i+0,j+1)b*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny;								
			//Direct
			h=8;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t2;
			//Inverse
			h=7;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t2;
	//THIRD NEIGHBOR
			/*Hoping (i,j)a<-->(i-2,j+0)*/						
			k1=((i-2+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			//Direct
			h=9;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t3;
			//Inverse
			h=10;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t3;
			/*Hoping (i,j)a<-->(i+2,j+0)a*/						
			k1=((i+2+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			//Direct
			h=10;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t3;
			//Inverse
			h=9;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t3;
			/*Hoping (i,j)a<-->(i+0,j-2)a*/						
			k1=((i+0+Nx)%Nx)*Ny+(j-2+Ny)%Ny;								
			//Direct
			h=11;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t3;
			//Inverse
			h=12;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t3;
			/*Hoping (i,j)a<-->(i+0,j+2)b*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+2+Ny)%Ny;								
			//Direct
			h=12;
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t3;
			//Inverse
			h=11;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )A					
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )B
			H.values[h+k1*v]		= t3;

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

		Floating pi2=((Floating)2)*((Floating)M_PI);
		cusp::complex<Floating> I(0.0f,1.0f);      
        //We determine the number of diagonals  v, which is related to the number of neighborns
        Integer h0,h1,k0,k1;      
		//Test if spin is used or not;
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
		//double phi=(((double)2)*((double)nPHI))/((double)Ny);
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){

			Scalar Phi[]	={cusp::exp(I*pi2*((Floating)(j*nPHI)/(Floating)Ny)),1};
//			std::cout<<Phi[0]<<" "<<Phi[1]<<std::endl;	
			//We start Filling the hoppings of the Matrix
			/*Hoping (i  ,j)--->(i+1,j)*/						
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;				
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			h0=1;														
			H.values[h0+k0*v]		= H.values[h0+k0*v]*Phi[0];
			/*Hoping (i+1,j)--->(i  ,j)*/
			h1=2;
			H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
			/*Hoping (i,j   )-->(i  ,j+1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
			h0=3;
			H.values[h0+k0*v]		= H.values[h0+k0*v]*Phi[1];
			/*Hoping (i ,j+1)-->(i,j)A*/
			h1=4;
			H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
		}
			
	}
	
template <typename Dim,typename Matrix>	void velocityx(Dim Nx, Dim Ny,Matrix& H,Matrix& V)
	{
	
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);      
		//FIRST NEIGHBORHS
		//           V2=(i,j+1)
		//             |
		//  V0=(i,j)_A |___ V3=(i+1,j )
		//             
		//               h=2
		//               |
		//          h=0  |___h=1

				//(i,j)-->(i+1,j  )B  ==>  (x,y)--> (x-1,y  )
				//(i,j)-->(i  ,j+1)B  ==>  (x,y)--> (x  ,y-1)

        Integer h0,h1,k0,k1;      
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
   		Scalar dr[] = {0,-1,1,0,0,-1,1,-1,1,-2,2,0,0};
		/**********FIRST NEIGHBORNS******************/
		V=H; 		
		for(int i=0;i<H.num_entries;i++) V.values[i]=0;
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
	//FIRST NEIGHBOR
			/*Hoping (i,j)a<-->(i+0,j+0)b*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			//Direct
			h0=1;
			V.values[h0+k0*v]		= I*dr[h0]*H.values[h0+k0*v];
			//Inverse
			V.values[h1+k1*v]		= cusp::conj(V.values[h0+k0*v]);;
			/*Hoping (i,j)a<-->(i+1,j+0)b*/						
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			//Direct
			h0=2;
			V.values[h0+k0*v]		= I*dr[h0]*H.values[h0+k0*v];
			//Inverse
			V.values[h1+k1*v]		= cusp::conj(V.values[h0+k0*v]);;
			/*Hoping (i,j)a<-->(i+0,j+1)b*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			//Direct
			h0=3;
			V.values[h0+k0*v]		= I*dr[h0]*H.values[h0+k0*v];
			//Inverse
			V.values[h1+k1*v]		= cusp::conj(V.values[h0+k0*v]);;
			/*Hoping (i,j)a<-->(i+1,j+1)b*/						
			k1=((i+1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			//Direct
			h0=4;
			V.values[h0+k0*v]		= I*dr[h0]*H.values[h0+k0*v];
			//Inverse
			V.values[h1+k1*v]		= cusp::conj(V.values[h0+k0*v]);;
	//SECOND NEIGHBOR
			/*Hoping (i,j)a<-->(i-1,j+0)*/						
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			//Direct
			h=5;
			H.values[h+k0*v]		= t2;
			//Inverse
			h=6;
			H.values[h+k1*v]		= t2;
			/*Hoping (i,j)a<-->(i+1,j+0)a*/						
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			//Direct
			h=6;
			H.values[h+k0*v]		= t2;
			//Inverse
			h=5;
			H.values[h+k1*v]		= t2;
			/*Hoping (i,j)a<-->(i+0,j-1)a*/						
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
			//Direct
			h=7;
			H.values[h+k0*v]		= t2;
			//Inverse
			h=8;
			H.values[h+k1*v]		= t2;
			/*Hoping (i,j)a<-->(i+0,j+1)b*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny;								
			//Direct
			h=8;
			H.values[h+k0*v]		= t2;
			//Inverse
			h=7;
			H.values[h+k1*v]		= t2;
	//THIRD NEIGHBOR
			/*Hoping (i,j)a<-->(i-2,j+0)*/						
			k1=((i-2+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			//Direct
			h=9;
			H.values[h+k0*v]		= t3;
			//Inverse
			h=10;
			H.values[h+k1*v]		= t3;
			/*Hoping (i,j)a<-->(i+2,j+0)a*/						
			k1=((i+2+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			//Direct
			h=10;
			H.values[h+k0*v]		= t3;
			//Inverse
			h=9;
			H.values[h+k1*v]		= t3;
			/*Hoping (i,j)a<-->(i+0,j-2)a*/						
			k1=((i+0+Nx)%Nx)*Ny+(j-2+Ny)%Ny;								
			//Direct
			h=11;
			H.values[h+k0*v]		= t3;
			//Inverse
			h=12;
			H.values[h+k1*v]		= t3;
			/*Hoping (i,j)a<-->(i+0,j+2)b*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+2+Ny)%Ny;								
			//Direct
			h=12;
			H.values[h+k0*v]		= t3;
			//Inverse
			h=11;
			H.values[h+k1*v]		= t3;
		}
		}
	
template <typename Dim,typename Matrix>	void velocityy(Dim Nx, Dim Ny,Matrix& H,Matrix& V)
	{
	
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;							//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
		cusp::complex<Floating> I(0.0f,1.0f);      
		//FIRST NEIGHBORHS
		//           V2=(i,j+1)
		//             |
		//  V0=(i,j)_A |___ V3=(i+1,j )
		//             
		//               h=2
		//               |		
		//          h=0  |___h=1

				//(i,j)-->(i+1,j  )B  ==>  (x,y)--> (x-1,y  )
				//(i,j)-->(i  ,j+1)B  ==>  (x,y)--> (x  ,y-1)

        Integer h0,h1,k0,k1;      
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
		/**********FIRST NEIGHBORNS******************/
		Scalar dr[]={ 0, 1};
		V=H; 		
		for(int i=0;i<H.num_entries;i++) V.values[i]=0;
	
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
			/*Hoping (i,j) -->(i-1,j  ) */						
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			//Direct
			h0=1;
			V.values[h0+k0*v]		= I*dr[0]*H.values[h0+k0*v];
			//Reversed
			h1=2;
			V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j) -->(i  ,j-1) */
			k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny;								
			//Direct
			h0=3;
			V.values[h0+k0*v]		 = I*dr[1]*H.values[h0+k0*v];
			//Reversed
			h1=4;
			V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
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
		}

	}
    		
    //End os the space square::
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

