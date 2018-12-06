
#ifndef MAT_OPERATOR_H
#define MAT_OPERATOR_H

#include "lattice_geometry.h"
#include "hamiltonian_parameters.h"


void FillHamiltonian(const int dimk, Complex* H,const HamParams& ham_params, const Real* k)
{
	const Complex E0[2]   ={ ham_params.E0[0] , ham_params.E0[1] };
	const Complex t1      = ham_params.t1;
	const Complex t2      = ham_params.t2 ;
	const Complex VI[2]   ={ ham_params.VI[0], ham_params.VI[1] };
	const Complex VR      =  ham_params.VR ;
	const Complex VPIA[2] = {ham_params.VPIA[0],ham_params.VPIA[1]};

	for(int n=0; n < dimk*dimk ; n++ )
		H[n]=0;
	
	const int i0=0, i1=0;
	//Go through the spin and the lattice indexes
	for(int io=0; io < NORB ; io++ )
	for(int is=0; is < SPIN ; is++ )
	{
		//Compute the sign associate with each index
		const Real	
		ss= 1 - 2*is,
		oo= 1 - 2*io ;
		//compute the final orbital/spin
		const int
		jo=1-io,
		js=1-is;
		//define the indexes for the nearest neighbors
		int DI0nn[3]={0,oo,0 };
		int DI1nn[3]={0,0 ,oo};
		//define the indexes of the next nearest neighbors
		int DI0nnn[6]={-1,+1, 0, 0 ,-1 ,+1 };
		int DI1nnn[6]={ 0, 0,-1,+1 ,+1 ,-1 };
		//define the inital index
		const int
		ki= is*NORB + io;
		const Real 
		ri[3]={ i0*A[0][0] + i1*A[1][0],
				i0*A[0][1] + i1*A[1][1], 
				 0 
			    };			    
//----------------------------ONSITE ENERGY/----------------------------/
		{
			const int
			kj= is*NORB + io;
			H [ki*dimk+kj]= E0[io] ;		
		}
//----------------------------NEAREST NEIGHBORS/----------------------------/
		for(int n=0;n<3;n++)
		{	
			const int
			kj= is*NORB + jo;
			//compute the r_j position
			const Real 
			rj[2]={ (i0+DI0nn[n])*A[0][0] + (i1+DI1nn[n])*A[1][0],
					(i0+DI0nn[n])*A[0][1] + (i1+DI1nn[n])*A[1][1] 
					};	
			//The basis vector is already taken into account, that is 
			//why it is removed below by  summing it
			const Complex
			val=t1*exp( I *( (rj[0]-ri[0])*k[0] +  (rj[1]-ri[1])*k[1] ));
			H [ki*dimk+kj]+= val ;
		}	
//----------------------------NEXT NEAREST NEIGHBORS/----------------------------/			
		const Complex
		iso_val[6]={(t2+ I*ss*oo*VI[io]/3./sqrt(3)),
					(t2- I*ss*oo*VI[io]/3./sqrt(3)),
					(t2- I*ss*oo*VI[io]/3./sqrt(3)),
					(t2+ I*ss*oo*VI[io]/3./sqrt(3)),
					(t2- I*ss*oo*VI[io]/3./sqrt(3)),
					(t2+ I*ss*oo*VI[io]/3./sqrt(3))
				  };
		for(int n=0;n<6;n++)
		{	
			const int
			kj= is*NORB + io;
			const Real 
			rj[2]={ ( i0+DI0nnn[n] )*A[0][0] + ( i1+DI1nnn[n] )*A[1][0],
					( i0+DI0nnn[n] )*A[0][1] + ( i1+DI1nnn[n] )*A[1][1] 
					};							
			const Complex
			val=iso_val[n]*exp( I *( (rj[0]-ri[0])*k[0] +  (rj[1]-ri[1])*k[1] ) );
			H [ki*dimk+kj]+=val ;
		}
//---------------------------- NEAREST NEIGHBORS SPIN FLIP/----------------------------/
		for(int n=0;n<3;++n)
		{
			const int
			kj= js*NORB + jo;
			//compute the r_j position
			const Real 
			rj[2]={ (i0+DI0nn[n])*A[0][0] + (i1+DI1nn[n])*A[1][0],
					(i0+DI0nn[n])*A[0][1] + (i1+DI1nn[n])*A[1][1] 
					};	
			//compute the differen between ri and rj
			const Real 
			del_r[3]={ rj[0]-ri[0]+(jo-io)*Delta[0] ,rj[1]-ri[1] + (jo-io)*Delta[1],0.  };
			//compute the rashba value
			const Complex 
			PhiRSO= 2.*I*VR*CrossProductDotZ(del_r, is,js )/3.,
			val   = PhiRSO*exp( I *( (rj[0]-ri[0])*k[0] +  (rj[1]-ri[1])*k[1] ) );
			H [ki*dimk+kj]+= val ;
		}

//---------------------------- NEAREST NEIGHBORS SPIN FLIP/-----------------------
		//------ PIA INTERACTION------
		for(int n=0;n<6;++n)
		{
			const int
			kj= js*NORB + io;
			const Real 
			rj[2]={ ( i0+DI0nnn[n] )*A[0][0] + ( i1+DI1nnn[n] )*A[1][0],
					( i0+DI0nnn[n] )*A[0][1] + ( i1+DI1nnn[n] )*A[1][1] 
					};						
			//compute the differen between ri and rj					
			const Real 
			del_r[3]={ rj[0]-ri[0] ,rj[1]-ri[1] ,rj[2]-ri[2]  };
			//Compute Pia phase
			const Complex 
			PhiPIA= 2.*I*VPIA[io]*CrossProductDotZ(del_r, is,js )/3,
			val   =PhiPIA*exp( I *( (rj[0]-ri[0])*k[0] +  (rj[1]-ri[1])*k[1] ) );
			H [ki*dimk+kj]+= val ;
		}
	}	
	
};


#endif
