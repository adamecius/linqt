



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
		std::cout<<"Creating the lattice for "<<v<<" Neighbors\n";
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
		if(v>=4){
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
		}
		if(v>=10){
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
		if(v>=13){
	//THIRD NEIGHBOR
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			/*Hoping (i,j)A-->(i+1,j+1)B*/						
			k1=((i+1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=10;														
			//DIRECT
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A					
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j  )B
			H.values[h+k0*v]		= t3;
			//INVERSE
			h=10;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j  )B				
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )A
			H.values[h+k1*v]		= t3;
			/*Hoping (i,j)A-->(i-1,j+1)B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h=11;
			//DIRECT
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			H.column_indices[h+k0*v]= k1;								//Site: (i+1,j  )B
			H.values[h+k0*v]		= t3;
			//INVERSE
			h=11;
			H.row_indices[h+k1*v]	= k1;								//Site: (i+1,j  )B				
			H.column_indices[h+k1*v]= k0  ;								//Site: (i  ,j  )A
			H.values[h+k1*v]		= t3;
			/*Hoping (i,j)A-->(i+1,j-1)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h=12;
			//DIRECT
			H.row_indices[h+k0*v]	= k0;								//Site: (i  ,j  )A				
			H.column_indices[h+k0*v]= k1;								//Site: (i  ,j+1)B
			H.values[h+k0*v]		= t3;
			//INVERSE
			h=12;
			H.row_indices[h+k1*v]	= k1;								//Site: (i  ,j+1)B				
			H.column_indices[h+k1*v]= k0;								//Site: (i  ,j  )A
			H.values[h+k1*v]		= t3;
		}		
		}    
	}

	namespace impurities{

		template <typename Dim,typename Matrix,typename FloatType,typename Vector> 
		void HtypeDistribution(Dim Nx, Dim Ny,Matrix& H, const FloatType p, const typename Matrix::value_type& U, const Vector& V,const FloatType eps0)
		{
				
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;						//This one will be use to counting and dimension parameters
        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.

        Integer D,h0,h1,k0,k1;      
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
		D=Nx*Ny;
			/*As shown in the following draw
			*	    V1=(i,j+1)_B ____ V2=(i,j+1)_A
			*					/    \
			*		V0=(i,j)_A /      \V3=(i+1,j+1)_B
			*				   \      /
			*	    V5=(i+1,j)_B\____/V4=(i+1,j)_A
			*
			* the vector V[6]=(V0,V1,V2,V3,V4,V5) has as components the hybridization amplitude
			* of the atoms in the ring within the Anderson Impurity model, 
			* in order to calculate the effective hopping tau[i][j], we perfomed a Canonical Transformation
			* in order to obtain an graphene only effective hamiltonian, within this approximation we obtain
			* tau_ij= - V_i conj(V_j)/ eps_0 where eps_0 is the energy of the localized level in the impurity in this case 
			* we are absorbing this inside V_i
			*/ 
			Scalar T[6][6];
			for(int i=0;i<6;i++)
				for(int j=0;j<6;j++)
					T[i][j]=V[i]*cusp::conj(V[j])*eps0;

		//We start by filling the system with the impurities
		std::cout<<"Setting up the impurities with "<<p<<" Concentration\n";
		for(Integer i=0;i<Nx;i++)
			for(Integer j=0;j<Ny;j++){
				//We first decide if the impurity is placed or not
				// with a probability p
				if((Floating)rand()/(Floating)RAND_MAX<p){
					//if selected randomly select whats the orientation
					//of the orbital symmetry by choosing a number between 
					//1 and 6
					int i0=rand()%6;
					//We then fill the hoppings
				//ATOM 0:  (i+0,j+0)A
						{
						//ZERO NEIGHBOR
						/*(i+0,j+0)A-->(i+0,j+0)A*/
						k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;				
						h0=0;	
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+i0)%6][(0+i0)%6]+U;
						//FIRST NEIGHBOR
						/*(i+0,j+0)A <--> (i+1,j+0)B*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;				
							//DIRECT								
						h0=2;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+i0)%6][(5+i0)%6];
							//INVERSE
						h1=2;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						/*(i+0,j+0)A <--> (i+0,j+1)B*/
						k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;				
							//DIRECT
						h0=3;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+i0)%6][(1+i0)%6];
							//INVERSE
						h1=3;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						//SECOND NEIGHBOR
						/*(i+0,j+0)A <--> (i+1,j+0)A*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;				
							//DIRECT								
						h0=5;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+i0)%6][(4+i0)%6];
							//INVERSE
						h1=4;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						/*(i+0,j+0)A <--> (i+0,j+1)A*/
						k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny;				
							//DIRECT
						h0=7;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+i0)%6][(2+i0)%6];
							//INVERSE
						h1=6;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						//THIRD NEIGHBOR
						/*(i+0,j+0)A <--> (i+1,j+1)B*/
						k1=((i+1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;				
							//DIRECT								
						h0=10;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+i0)%6][(3+i0)%6];
							//INVERSE
						h1=10;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
					}
				//ATOM 1:  (i+0,j+1)B
						{
						//ZERO NEIGHBOR
						/*(i+0,j+1)B-->(i+0,j+1)B*/
						k0=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;				
						h0=0;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+i0)%6][(1+i0)%6]+U;
						//FIRST NEIGHBOR
						/*(i+0,j+1)B-->(i+0,j+0)A*/
						/***This one was already set on ATOM 0 section***/
						/*(i+0,j+1)B-->(i+0,j+1)A*/
						k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny;				
							//DIRECT								
						h0=1;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+i0)%6][(2+i0)%6];
							//INVERSE
						h1=1;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						//SECOND NEIGHBOR
						/*(i+0,j+1)B-->(i+1,j+1)B*/
						k1=((i+1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;				
							//DIRECT								
						h0=5;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+i0)%6][(3+i0)%6];
							//INVERSE
						h1=4;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						/*(i+0,j+1)B-->(i+1,j+0)B*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;				
							//DIRECT
						h0=9;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+i0)%6][(5+i0)%6];
							//INVERSE
						h1=8;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						//THIRD NEIGHBOR
						/*(i+0,j+1)B-->(i+1,j+0)A*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;				
							//DIRECT								
						h0=12;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+i0)%6][(4+i0)%6];
							//INVERSE
						h1=12;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
					}
				//ATOM 2:  (i+0,j+1)A
						{
						//ZERO NEIGHBOR
						/*(i+0,j+1)A-->(i+0,j+1)A*/
						k0=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny ;				
						h0=0;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(2+i0)%6][(2+i0)%6]+U;
						//FIRST NEIGHBOR
						/*(i+0,j+1)A-->(i+0,j+1)B*/
						/***This one was already set on ATOM 1 section***/
						/*(i+0,j+1)A-->(i+1,j+1)B*/
						k1=((i+1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;				
							//DIRECT
						h0=2;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(2+i0)%6][(3+i0)%6];
							//INVERSE
						h1=2;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						//SECOND NEIGHBOR
						/*(i+0,j+1)A-->(i+0,j+0)A*/
						/***This one was already set on ATOM 0 section***/
						/*(i+0,j+1)A-->(i+1,j+0)A*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;				
							//DIRECT
						h0=9;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(2+i0)%6][(4+i0)%6];
							//INVERSE
						h1=8;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						//THIRD NEIGHBOR
						/*(i+0,j+1)A-->(i+1,j+0)B*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;				
							//DIRECT								
						h0=12;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(2+i0)%6][(5+i0)%6];
							//INVERSE
						h1=12;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
					}
				//ATOM 3:  (i+1,j+1)B
						{
						//ZERO NEIGHBOR
						/*(i+1,j+1)B-->(i+1,j+1)B*/
						k0=((i+1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;				
						h0=0;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(3+i0)%6][(3+i0)%6]+U;
						//FIRST NEIGHBOR
						/*(i+1,j+1)B-->(i+0,j+1)A*/
						/***This one was already set on ATOM 2 section***/
						/*(i+1,j+1)B-->(i+1,j+0)A*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;				
							//DIRECT
						h0=3;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(3+i0)%6][(4+i0)%6];
							//INVERSE
						h1=3;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						//SECOND NEIGHBOR
						/*(i+1,j+1)B-->(i+0,j+1)B*/
						/***This one was already set on ATOM 1 section***/
						/*(i+1,j+1)B-->(i+1,j+0)B*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;				
							//DIRECT
						h0=6;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(3+i0)%6][(5+i0)%6];
							//INVERSE
						h1=7;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						//THIRD NEIGHBOR
						/*(i+1,j+1)B-->(i+0,j+0)a*/
						/***This one was already set on ATOM 0 section***/
					}
				//ATOM 4:  (i+1,j+0)A
						{
						//ZERO NEIGHBOR
						/*(i+1,j+0)A-->(i+1,j+0)A*/
						k0=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny;				
						h0=0;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(4+i0)%6][(4+i0)%6]+U;
						//FIRST NEIGHBOR
						/*(i+1,j+0)A-->(i+1,j+1)B*/
						/***This one was already set on ATOM 3 section***/
						/*(i+1,j+0)A-->(i+1,j+0)B*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;				
							//DIRECT
						h0=1;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(4+i0)%6][(5+i0)%6];
							//INVERSE
						h1=1;														
						H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
						//SECOND NEIGHBOR
						/*(i+1,j+0)A-->(i+0,j+0)A*/
						/***This one was already set on ATOM 0 section***/
						/*(i+1,j+0)A-->(i+0,j+1)A*/
						/***This one was already set on ATOM 2 section***/
						//THIRD NEIGHBOR
						/*(i+1,j+0)A-->(i+0,j+1)B*/
						/***This one was already set on ATOM 1 section***/
					}
				//ATOM 5:  (i+1,j+0)B
						{
						//ZERO NEIGHBOR
						/*(i+1,j+0)B-->(i+1,j+0)B*/
						k0=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;				
						h0=0;														
						H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(5+i0)%6][(5+i0)%6]+U;
						//FIRST NEIGHBOR
						/*(i+1,j+0)B-->(i+0,j+0)A*/
						/***This one was already set on ATOM 0 section***/
						/*(i+1,j+0)B-->(i+1,j+0)A*/
						/***This one was already set on ATOM 4 section***/
						//SECOND NEIGHBOR
						/*(i+1,j+0)B-->(i+0,j+1)B*/
						/***This one was already set on ATOM 1 section***/
						/*(i+1,j+0)B-->(i+1,j+1)B*/
						/***This one was already set on ATOM 3 section***/
						//THIRD NEIGHBOR
						/*(i+1,j+0)B-->(i+0,j+1)A*/
						/***This one was already set on ATOM 2 section***/
					}
				}
			}
		}

		template <typename Dim,typename Matrix,typename FloatType,typename Vector> 
		void Htype(Dim Nx, Dim Ny,Matrix& H, const int i0, const int j0, const typename Matrix::value_type& U, const Vector& V,const FloatType eps0)
		{				
        //We pass the variable type from the Matrix to the local parameters
			typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
	        typedef typename Matrix::index_type Integer;						//This one will be use to counting and dimension parameters
	        typedef typename Matrix::value_type::value_type Floating;			//This one will define if the code should run on the device or the host.
	        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
	        Integer D,h0,h1,k0,k1;      
	        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
			D=Nx*Ny;
	        /********************************Pristine Hamiltonian***********************************/
			/*Graphene's Lattices vector are necessary to implement this method
			 * because we need to measure the distance between atoms in order to 
			 * incorporate the constraing of the circle */
			 Floating a1[]={ 1.0*sqrt(3), 0.0};
			 Floating a2[]={ 0.5*sqrt(3), 1.5};
			 Floating d1[]={-0.5*sqrt(3),-0.5};
			/***Center of the Disk initially fixed at the center of graphene. CAN BE EXTENDED***/
			Floating r []={a1[0]*i0+a2[0]*j0,a1[1]*i0+a2[1]*j0};
			std::cout<<i0<<" "<<j0<<std::endl;
			std::cout<<r[0]<<" "<<r[1]<<std::endl;
			/*As shown in the following draw
			*	    V1=(i,j+1)_B ____ V2=(i,j+1)_A
			*					/    \
			*		V0=(i,j)_A /      \V3=(i+1,j+1)_B
			*				   \      /
			*	    V5=(i+1,j)_B\____/V4=(i+1,j)_A
			*
			* the vector V[6]=(V0,V1,V2,V3,V4,V5) has as components the hybridization amplitude
			* of the atoms in the ring within the Anderson Impurity model, 
			* in order to calculate the effective hopping tau[i][j], we perfomed a Canonical Transformation
			* in order to obtain an graphene only effective hamiltonian, within this approximation we obtain
			* tau_ij= - V_i conj(V_j)/ eps_0 where eps_0 is the energy of the localized level in the impurity in this case 
			* we are absorbing this inside V_i
			*/ 
			Scalar T[6][6];
			for(int m=0;m<6;m++)
				for(int n=0;n<6;n++){
					T[m][n]=V[m]*cusp::conj(V[n])*eps0;
					std::cout<<T[m][n]<<" "<<std::endl;
					}
			//if selected randomly select whats the orientation
			//of the orbital symmetry by choosing a number between 
			//1 and 6
			int l0=rand()%6;
			//We then fill the hoppings
			//ATOM 0:  (i+0,j+0)A
			{
				//ZERO NEIGHBOR
				/*(i+0,j+0)A-->(i+0,j+0)A*/
				k0=((i0+0+Nx)%Nx)*Ny+(j0+0+Ny)%Ny;				
				h0=0;	
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+l0)%6][(0+l0)%6]+U;
				//FIRST NEIGHBOR
				/*(i+0,j+0)A <--> (i+1,j+0)B*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny+D;
				//DIRECT
				h0=2;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+l0)%6][(5+l0)%6];
				//INVERSE
				h1=2;
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				/*(i+0,j+0)A <--> (i+0,j+1)B*/
				k1=((i0+0+Nx)%Nx)*Ny+(j0+1+Ny)%Ny+D;
				//DIRECT
				h0=3;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+l0)%6][(1+l0)%6];
				//INVERSE
				h1=3;
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				//SECOND NEIGHBOR
				/*(i+0,j+0)A <--> (i+1,j+0)A*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny;
				//DIRECT
				h0=5;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+l0)%6][(4+l0)%6];
				//INVERSE
				h1=4;
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				/*(i+0,j+0)A <--> (i+0,j+1)A*/
				k1=((i0+0+Nx)%Nx)*Ny+(j0+1+Ny)%Ny;
				//DIRECT
				h0=7;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+l0)%6][(2+l0)%6];
				//INVERSE
				h1=6;
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				//THIRD NEIGHBOR
				/*(i+0,j+0)A <--> (i+1,j+1)B*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+1+Ny)%Ny+D;
				//DIRECT
				h0=10;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(0+l0)%6][(3+l0)%6];
				//INVERSE
				h1=10;
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
			}
		//ATOM 1:  (i+0,j+1)B
			{
				//ZERO NEIGHBOR
				/*(i+0,j+1)B-->(i+0,j+1)B*/
				k0=((i0+0+Nx)%Nx)*Ny+(j0+1+Ny)%Ny+D;
				h0=0;	
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+l0)%6][(1+l0)%6]+U;
				//FIRST NEIGHBOR
				/*(i+0,j+1)B-->(i+0,j+0)A*/
				/***This one was already set on ATOM 0 section***/
				/*(i+0,j+1)B-->(i+0,j+1)A*/
				k1=((i0+0+Nx)%Nx)*Ny+(j0+1+Ny)%Ny;
				//DIRECT
				h0=1;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+l0)%6][(2+l0)%6];
				//INVERSE
				h1=1;
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				//SECOND NEIGHBOR
				/*(i+0,j+1)B-->(i+1,j+1)B*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+1+Ny)%Ny+D;
				//DIRECT
				h0=5;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+l0)%6][(3+l0)%6];
				//INVERSE
				h1=4;
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				/*(i+0,j+1)B-->(i+1,j+0)B*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny+D;
				//DIRECT
				h0=9;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+l0)%6][(5+l0)%6];
				//INVERSE
				h1=8;
				//H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				//THIRD NEIGHBOR
				/*(i+0,j+1)B-->(i+1,j+0)A*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny;
				//DIRECT
				h0=12;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(1+l0)%6][(4+l0)%6];
				//INVERSE
				h1=12;
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
			}
		//ATOM 2:  (i+0,j+1)A
			{
				//ZERO NEIGHBOR
				/*(i+0,j+1)A-->(i+0,j+1)A*/
				k0=((i0+0+Nx)%Nx)*Ny+(j0+1+Ny)%Ny ;
				h0=0;
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(2+l0)%6][(2+l0)%6]+U;
				//FIRST NEIGHBOR
				/*(i+0,j+1)A-->(i+0,j+1)B*/
				/***This one was already set on ATOM 1 section***/
				/*(i+0,j+1)A-->(i+1,j+1)B*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+1+Ny)%Ny+D;				
				//DIRECT
				h0=2;														
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(2+l0)%6][(3+l0)%6];
				//INVERSE
				h1=2;														
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				//SECOND NEIGHBOR
				/*(i+0,j+1)A-->(i+0,j+0)A*/
				/***This one was already set on ATOM 0 section***/
				/*(i+0,j+1)A-->(i+1,j+0)A*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny;				
				//DIRECT
				h0=9;														
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(2+l0)%6][(4+l0)%6];
				//INVERSE
				h1=8;														
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				//THIRD NEIGHBOR
				/*(i+0,j+1)A-->(i+1,j+0)B*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny+D;				
				//DIRECT								
				h0=12;														
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(2+l0)%6][(5+l0)%6];
				//INVERSE
				h1=12;														
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
			}
		//ATOM 3:  (i+1,j+1)B
			{
				//ZERO NEIGHBOR
				/*(i+1,j+1)B-->(i+1,j+1)B*/
				k0=((i0+1+Nx)%Nx)*Ny+(j0+1+Ny)%Ny+D;				
				h0=0;														
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(3+l0)%6][(3+l0)%6]+U;
				//FIRST NEIGHBOR
				/*(i+1,j+1)B-->(i+0,j+1)A*/
				/***This one was already set on ATOM 2 section***/
				/*(i+1,j+1)B-->(i+1,j+0)A*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny;				
				//DIRECT
				h0=3;														
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(3+l0)%6][(4+l0)%6];
				//INVERSE
				h1=3;														
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				//SECOND NEIGHBOR
				/*(i+1,j+1)B-->(i+0,j+1)B*/
				/***This one was already set on ATOM 1 section***/
				/*(i+1,j+1)B-->(i+1,j+0)B*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny+D;				
					//DIRECT
				h0=6;														
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(3+l0)%6][(5+l0)%6];
					//INVERSE
				h1=7;														
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				//THIRD NEIGHBOR
				/*(i+1,j+1)B-->(i+0,j+0)a*/
				/***This one was already set on ATOM 0 section***/
			}
		//ATOM 4:  (i+1,j+0)A
			{
				//ZERO NEIGHBOR
				/*(i+1,j+0)A-->(i+1,j+0)A*/
				k0=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny;				
				h0=0;														
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(4+l0)%6][(4+l0)%6]+U;
				//FIRST NEIGHBOR
				/*(i+1,j+0)A-->(i+1,j+1)B*/
				/***This one was already set on ATOM 3 section***/
				/*(i+1,j+0)A-->(i+1,j+0)B*/
				k1=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny+D;				
					//DIRECT
				h0=1;														
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(4+l0)%6][(5+l0)%6];
					//INVERSE
				h1=1;														
				H.values[h1+k1*v]		= cusp::conj(H.values[h0+k0*v]);
				//SECOND NEIGHBOR
				/*(i+1,j+0)A-->(i+0,j+0)A*/
				/***This one was already set on ATOM 0 section***/
				/*(i+1,j+0)A-->(i+0,j+1)A*/
				/***This one was already set on ATOM 2 section***/
				//THIRD NEIGHBOR
				/*(i+1,j+0)A-->(i+0,j+1)B*/
				/***This one was already set on ATOM 1 section***/
			}
		//ATOM 5:  (i+1,j+0)B
			{
				//ZERO NEIGHBOR
				/*(i+1,j+0)B-->(i+1,j+0)B*/
				k0=((i0+1+Nx)%Nx)*Ny+(j0+0+Ny)%Ny+D;				
				h0=0;														
				H.values[h0+k0*v]		= H.values[h0+k0*v]+T[(5+l0)%6][(5+l0)%6]+U;
				//FIRST NEIGHBOR
				/*(i+1,j+0)B-->(i+0,j+0)A*/
				/***This one was already set on ATOM 0 section***/
				/*(i+1,j+0)B-->(i+1,j+0)A*/
				/***This one was already set on ATOM 4 section***/
				//SECOND NEIGHBOR
				/*(i+1,j+0)B-->(i+0,j+1)B*/
				/***This one was already set on ATOM 1 section***/
				/*(i+1,j+0)B-->(i+1,j+1)B*/
				/***This one was already set on ATOM 3 section***/
				//THIRD NEIGHBOR
				/*(i+1,j+0)B-->(i+0,j+1)A*/
				/***This one was already set on ATOM 2 section***/
			}
		
		}
	
		template <typename Dim, typename Matrix, typename FloatType>
		void TtypeDistribution(Dim Nx, Dim Ny,Matrix& H, const FloatType p,const typename Matrix::value_type& U,const typename Matrix::value_type& tU)
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


	}	

	namespace potentials{
		template <typename Dim, typename Matrix, typename Floating>
		void CenteredDoubleDisk(Dim Nx, Dim Ny,  Matrix& H, const Floating Umin, const Floating Rmin, const Floating Umax, const Floating Rmax)
		{
        //We pass the variable type from the Matrix to the local parameters
        typedef typename Matrix::value_type Scalar;							//This one will be use to the scalar parameters (i.e vector and matrix componnents
        typedef typename Matrix::index_type Integer;						//This one will be use to counting and dimension parameters
        typedef typename Matrix::memory_space MemorySpace;					//This one will define if the code should run on the device or the host.
        Integer i0,j0,k0,h0,D;
        Floating radius;

		//This method is design to incorpored a double disk inside graphene
		D=Nx*Ny;
        Integer v=(Integer)((double)H.num_entries/(double)H.num_cols);
        /********************************Pristine Hamiltonian***********************************/
		/*Graphene's Lattices vector are necessary to implement this method
		 * because we need to measure the distance between atoms in order to 
		 * incorporate the constraing of the circle */
		 Floating a1[]={ 1.0*sqrt(3), 0.0};
		 Floating a2[]={ 0.5*sqrt(3), 1.5};
		 Floating d1[]={-0.5*sqrt(3),-0.5};
		/***Center of the Disk initially fixed at the center of graphene. CAN BE EXTENDED***/
		i0=(Integer)((Floating)Nx/2.0); /*Centered at the X direction*/
		j0=(Integer)((Floating)Ny/2.0); /*Centered at the Y direction*/
		Floating rc[]={a1[0]*i0+a2[0]*j0,a1[1]*i0+a2[1]*j0};
		Floating r []={a1[0]*i0+a2[0]*j0,a1[1]*i0+a2[1]*j0};
		/***Disk Parameters. CAN BE EXTENDED***/
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++)for(Integer z=0;z<2;z++){
            /***********************RED A------->X*****************************/
			//*Calculate the r point associated with (i,j)_z indexes
			r[0]=a1[0]*i+a2[0]*j+d1[0]*z-rc[0];
			r[1]=a1[1]*i+a2[1]*j+d1[1]*z-rc[1];
			/*Check if the point belongs to the selected region
			 * In this case the intersecption between the two circles C1=(r0,Rmax) and C2=(r0,Rmin)
			 */
			radius=sqrt(r[0]*r[0]+r[1]*r[1]);
			if(radius < Rmin){
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+z*D;
				h0=0;
				H.values[h0+k0*v]		=H.values[h0+k0*v]+Umin;
				}else 
				if( radius <Rmax){
					k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+z*D;
					h0=0;
					H.values[h0+k0*v]	=H.values[h0+k0*v]+Umax;;
					}
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

				//Minimal set of hoppings
				//h=1;	(i,j)A-->(i+0,j+0)B  ==>  (x,y)--> (x-sqrt(3)/2	,y-1/2	)
				//h=2;	(i,j)A-->(i+1,j+0)B  ==>  (x,y)--> (x+sqrt(3)/2	,y-1/2	)
				//h=3;	(i,j)A-->(i+0,j+1)B  ==>  (x,y)--> (x+0			,y+1	)
				//h=4;	(i,j)X-->(i-1,j  )X  ==>  (x,y)--> (x-sqrt(3)	,y+0	)
				//h=6;	(i,j)X-->(i  ,j-1)X  ==>  (x,y)--> (x-sqrt(3)/2	,y-3/2	)
				//h=8;	(i,j)X-->(i-1,j+1)X  ==>  (x,y)--> (x-sqrt(3)/2	,y+3/2	)
				//h=10;	(i,j)A-->(i+1,j+1)B  ==>  (x,y)--> (x+sqrt(3)	,y+1	)
				//h=11;	(i,j)A-->(i-1,j+1)B  ==>  (x,y)--> (x-sqrt(3)	,y+1	)
				//h=11;	(i,j)A-->(i+1,j-1)B  ==>  (x,y)--> (x+0			,y-3/2	)

        Integer D,h0, h1,k0,k1;      
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
        D=Nx*Ny;
		/**********FIRST NEIGHBORNS******************/
		Scalar 		dr[]={	/**NN  **/-0.5*sqrt(3), 0.5*sqrt(3), 0.0			,
							/**NNN **/-1.0*sqrt(3),-0.5*sqrt(3),-0.5*sqrt(3)	,
							/**NNNN**/ 1.0*sqrt(3),-1.0*sqrt(3), 0.0			};
		V=H; for(int i=0;i<H.num_entries;i++) V.values[i]=0;
		
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
		if(v>=3){
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;
			//FIRST NEIGHBORS
			/*Hoping (i,j)A-->(i,j)B*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			//Direct
			h0=1;
			V.values[h0+k0*v]		 = I*dr[0]*H.values[h0+k0*v];
			/*Hoping (i,j)B-->(i,j)A*/
			//Inverse
			h1=1;
			V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i+1,j)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			//Direct
			h0=2;														
			V.values[h0+k0*v]		 = I*dr[1]*H.values[h0+k0*v];
			//Inverse
			h1=2;
			V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i,j+1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			//Direct
			h0=3;
			V.values[h0+k0*v]		 = I*dr[2]*H.values[h0+k0*v];
			//Inverse
			h1=3;
			V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
			}
		if(v>=10){
			for(int z0=0;z0<2;z0++){
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+z0*D;	
				//Second Neighborns A Lattice->Lattice
				/*Hoping (i,j)A-->(i-1,j)A*/
				k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+z0*D;								
				//Direct
				h0=4;														
				V.values[h0+k0*v]		 = I*dr[3]*H.values[h0+k0*v];
				//Inverse
				h1=5;
				V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i,j-1)A*/
				k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny+z0*D;								
				//Direct
				h0=6;														
				V.values[h0+k0*v]		 = I*dr[4]*H.values[h0+k0*v];
				//Inverse
				h1=7;
				V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
				/*Hoping (i,j)A-->(i,j+1)B*/
				k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+z0*D;								
				//Direct
				h0=8;														
				V.values[h0+k0*v]		 = I*dr[5]*H.values[h0+k0*v];
				//Inverse
				h1=9;
				V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
				}
			}
		if(v>=13){
	//THIRD NEIGHBOR
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			/*Hoping (i,j)A-->(i+1,j+1)B*/						
			k1=((i+1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h0=10;														
			//DIRECT
			H.values[h0+k0*v]		= I*dr[6]*H.values[h0+k0*v];
			//INVERSE
			h1=10;
			H.values[h1+k1*v]		= cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i-1,j+1)B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h0=11;
			//DIRECT
			H.values[h0+k0*v]		= I*dr[7]*H.values[h0+k0*v];
			//INVERSE
			h1=11;
			H.values[h1+k1*v]		= cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i+1,j-1)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h0=12;
			//DIRECT
			H.values[h0+k0*v]		= I*dr[8]*H.values[h0+k0*v];
			//INVERSE
			h1=12;
			H.values[h1+k1*v]		= cusp::conj(V.values[h1+k0*v]);				
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

				//h=1;	(i,j)A-->(i+0,j+0)B  ==>  (x,y)--> (x-sqrt(3)/2	,y-1/2	)
				//h=2;	(i,j)A-->(i+1,j+0)B  ==>  (x,y)--> (x+sqrt(3)/2	,y-1/2	)
				//h=3;	(i,j)A-->(i+0,j+1)B  ==>  (x,y)--> (x+0			,y+1	)
				//h=4;	(i,j)X-->(i-1,j  )X  ==>  (x,y)--> (x-sqrt(3)	,y+0	)
				//h=6;	(i,j)X-->(i  ,j-1)X  ==>  (x,y)--> (x-sqrt(3)/2	,y-3/2	)
				//h=8;	(i,j)X-->(i-1,j+1)X  ==>  (x,y)--> (x-sqrt(3)/2	,y+3/2	)
				//h=10;	(i,j)A-->(i+1,j+1)B  ==>  (x,y)--> (x+sqrt(3)	,y+1	)
				//h=11;	(i,j)A-->(i-1,j+1)B  ==>  (x,y)--> (x-sqrt(3)	,y+1	)
				//h=11;	(i,j)A-->(i+1,j-1)B  ==>  (x,y)--> (x+0			,y-3/2	)
        Integer D,h0, h1,k0,k1;      
        Integer v=(Integer)((double)H.num_entries/((double)H.num_cols));
        D=Nx*Ny;
		/**********FIRST NEIGHBORNS******************/
		Scalar dr[]		={	/**NN  **/-0.5,-0.5	, 1.0,
							/**NNN **/ 0.0,-1.5	, 1.5,
							/**NNNN**/ 1.0, 1.0	,-1.5};

		V=H; for(int i=0;i<H.num_entries;i++) V.values[i]=0;
		
		for(Integer i=0;i<Nx;i++)for(Integer j=0;j<Ny;j++){
			//We start Filling the hoppings of the Matrix
		if(v>=3){
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;
			//FIRST NEIGHBORS
			/*Hoping (i,j)A-->(i,j)B*/						
			k1=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			//Direct
			h0=1;
			V.values[h0+k0*v]		 = I*dr[0]*H.values[h0+k0*v];
			/*Hoping (i,j)B-->(i,j)A*/
			//Inverse
			h1=1;
			V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i+1,j)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+D;								
			//Direct
			h0=2;														
			V.values[h0+k0*v]		 = I*dr[1]*H.values[h0+k0*v];
			//Inverse
			h1=2;
			V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i,j+1)B*/
			k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			//Direct
			h0=3;
			V.values[h0+k0*v]		 = I*dr[2]*H.values[h0+k0*v];
			//Inverse
			h1=3;
			V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
			}
		if(v>=10){
			for(int z0=0;z0<2;z0++){
				k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+z0*D;	
				//Second Neighborns A Lattice->Lattice
				/*Hoping (i,j)A-->(i-1,j)A*/
				k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+z0*D;								
				//Direct
				h0=4;														
				V.values[h0+k0*v]		 = I*dr[3]*H.values[h0+k0*v];
				//Inverse
				h1=5;
				V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i,j-1)A*/
				k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny+z0*D;								
				//Direct
				h0=6;														
				V.values[h0+k0*v]		 = I*dr[4]*H.values[h0+k0*v];
				//Inverse
				h1=7;
				V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
				/*Hoping (i,j)A-->(i,j+1)B*/
				k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+z0*D;								
				//Direct
				h0=8;														
				V.values[h0+k0*v]		 = I*dr[5]*H.values[h0+k0*v];
				//Inverse
				h1=9;
				V.values[h1+k1*v]		 = cusp::conj(V.values[h0+k0*v]);
				}
			}
		if(v>=13){
	//THIRD NEIGHBOR
			k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny;								
			/*Hoping (i,j)A-->(i+1,j+1)B*/						
			k1=((i+1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h0=10;														
			//DIRECT
			H.values[h0+k0*v]		= I*dr[6]*H.values[h0+k0*v];
			//INVERSE
			h1=10;
			H.values[h1+k1*v]		= cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i-1,j+1)B*/
			k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+D;								
			h0=11;
			//DIRECT
			H.values[h0+k0*v]		= I*dr[7]*H.values[h0+k0*v];
			//INVERSE
			h1=11;
			H.values[h1+k1*v]		= cusp::conj(V.values[h0+k0*v]);
			/*Hoping (i,j)A-->(i+1,j-1)B*/
			k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny+D;								
			h0=12;
			//DIRECT
			H.values[h0+k0*v]		= I*dr[8]*H.values[h0+k0*v];
			//INVERSE
			h1=12;
			H.values[h1+k1*v]		= cusp::conj(V.values[h1+k0*v]);				
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
	
    template <typename Dim, typename Matrix>	void ISOa(Dim Nx, Dim Ny, Matrix& H,double lambda,double p)
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
        Floating ISOa=((Floating)1)*lambda/((Floating)3*sqrt(3));
		for(Integer i=0;i<Nx;i++)
			for(Integer j=0;j<Ny;j++)
				for(int z0=0;z0<2;z0++){
					int z1=1-z0;
					int zz=1-2*z0;					
					if((Floating)rand()/(Floating)RAND_MAX<p){
						k0=((i+0+Nx)%Nx)*Ny+(j+0+Ny)%Ny+z0*D;								
						/*Hoping (i,j)sz<-->(i-1,j  )sz*/
						k1=((i-1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+z0*D;								
						h0=4;
						H.values[h0+k0*v]		= H.values[h0+k0*v]-zz*I*ISOa;
						h1=5;
						H.values[h1+k1*v]		= H.values[h1+k1*v]+zz*I*ISOa;
						/*Second (i,j)sz<-->(i+1,j  )sz*/
						k1=((i+1+Nx)%Nx)*Ny+(j+0+Ny)%Ny+z0*D;								
						h0=5;
						H.values[h0+k0*v]		= H.values[h0+k0*v]+zz*I*ISOa;
						h1=4;
						H.values[h1+k1*v]		= H.values[h1+k1*v]-zz*I*ISOa;
						/*Second (i,j)sz<-->(i  ,j-1)sz*/
						k1=((i+0+Nx)%Nx)*Ny+(j-1+Ny)%Ny+z0*D;								
						h0=6;
						H.values[h0+k0*v]		= H.values[h0+k0*v]+zz*I*ISOa;
						h1=7;
						H.values[h1+k1*v]		= H.values[h1+k1*v]-zz*I*ISOa;
						/*Second (i,j)sz<-->(i  ,j+1)sz*/
						k1=((i+0+Nx)%Nx)*Ny+(j+1+Ny)%Ny+z0*D;								
						h0=7;
						H.values[h0+k0*v]		= H.values[h0+k0*v]-zz*I*ISOa;
						h1=6;
						H.values[h1+k1*v]		= H.values[h1+k1*v]+zz*I*ISOa;
						/*Second (i,j)sz<-->(i-1,j+1)sz*/
						k1=((i-1+Nx)%Nx)*Ny+(j+1+Ny)%Ny+z0*D;								
						h0=8;
						H.values[h0+k0*v]		= H.values[h0+k0*v]+zz*I*ISOa;
						h1=9;
						H.values[h1+k1*v]		= H.values[h1+k1*v]-zz*I*ISOa;
						/*Second (i,j)sz<-->(i+1,j-1)sz*/
						k1=((i+1+Nx)%Nx)*Ny+(j-1+Ny)%Ny+z0*D;								
						h0=9;
						H.values[h0+k0*v]		= H.values[h0+k0*v]-zz*I*ISOa;
						h1=8;
						H.values[h1+k1*v]		= H.values[h1+k1*v]+zz*I*ISOa;
						}
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
    std::cout<<"Refining the Matrix: Use of memory without refining "<<sizeof(H.values[0])*DD*v*1.0/1000000000<<"GB";
    std::cout<<" and after "<<sizeof(H.values[0])*nnz*1.0/1000000000<<"GB";
	std::cout<<" the memory economy was "<<100*(1.0-(double)nnz/(double)(DD*v))<<"%"<<std::endl;
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

