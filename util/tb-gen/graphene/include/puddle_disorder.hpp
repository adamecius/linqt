#ifndef PUDDLE_DISORDER
#define PUDDLE_DISORDER

#include<sys/time.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <sys/time.h>
#include <unistd.h>
#include <limits>

#include "lattice_geometry.h"


namespace puddle
{
int
IndexesToIndex(	const int i0, const int dim0,
		const int i1, const int dim1,
		const int i2, const int dim2,
		const int i3, const int dim3)
{
	return ( ( (i0+dim0)%dim0*dim1 + (i1+dim1)%dim1 )*dim2 + (i2+dim2)%dim2 )*dim3 + i3 ;
}
void
IndexToIndexes(	const int& k0,
				 int& i0, const int dim0,
				 int& i1, const int dim1,
				 int& i2, const int dim2,
				 int& i3, const int dim3)
{
	i3 =  k0%dim3;
	i2 = (k0/dim3)%dim2;
	i1 = (k0/dim3/dim2)%dim1;
	i0 = (k0/dim3/dim2/dim1)%dim0;
}
}				

int PuddleDisorder(	
		const int n0, const int n1, 
		const int norb, const int nspin,
		Real puddHeight,
		Real puddConcen,
		Real puddRange,	
		Real Ex,	
		std::vector<Real>& diagElements, Real seed)
{
	std::cout<<"I am using "<<(int)seed<<std::endl;
	srand((int)seed);
	const int dim=n0*n1*norb*nspin;
	const double epsilon=std::numeric_limits<double>::epsilon();


	//Check if the height or the concentration are not zero
	if(puddHeight <= epsilon || puddConcen  <= epsilon)
		return 0;

	///Estimate the cutoff indexes		
	const Real 
	L0=sqrt(pow(A[0][0],2.) + pow(A[0][1],2.) ) ,
	L1=sqrt(pow(A[1][0],2.) + pow(A[1][1],2.) ) ;

	//Evaluate how many unit cells I am going to evaluate the puddle function
	//it is bassically defined in terms of the puddle range
	const Real 
	CellCutOff = 20.*puddRange ; 

	const int 
	i0c=CellCutOff/L0  ,  //Dividing by the max distance gives you the 
	i1c=CellCutOff/L1  ;  // the number of cells
	std::cout<<"The puddle function are evaluated within the following ";
	std::cout<<"super cells : "<<i0c<<" "<<i1c<<std::endl;

	int imp_num=0;
	// first we go through all the possible positions
	for(int ip0=0; ip0 < n0 ; ip0++ )
	for(int ip1=0; ip1 < n1  ; ip1++ )
	for(int ipo=0; ipo < norb ; ipo++ )
	if(  (double)rand()/(double)RAND_MAX < puddConcen ) //Select a puddle center
	{
		const double 
		Height = puddHeight*(2.*( (float)rand()/(float)RAND_MAX ) - 1); ;

		//Check the neighboring unit cells
		for(int i0=-i0c; i0 <= i0c ; i0++ )
		for(int i1=-i1c; i1 <= i1c ; i1++ )
		for(int IO=0; IO < norb ; IO++ )
		{
			const int 
			I0 = ip0+i0,
			I1 = ip1+i1;

			const double
			rdiff[2] ={ (ip0-I0)*A[0][0] + (ip1-I1)*A[1][0]+ (IO-ipo)*Delta[0],
						(ip0-I0)*A[0][1] + (ip1-I1)*A[1][1]+ (IO-ipo)*Delta[1]
						},
			dist = 	sqrt( rdiff[0]*rdiff[0] + rdiff[1]*rdiff[1] );
				
			if( dist < CellCutOff*0.5 )
			{
				
				for( int IS=0; IS<nspin; IS++)
				{
					const double ss= 1.0 - 2.0*IS ;
					const double
					Ei= ( Height + ss*Ex)*exp( -0.5*dist*dist/puddRange/puddRange );
					const int
					k0= puddle::IndexesToIndex(	I0, n0,
												I1, n1,
												IO, norb,
												IS, nspin );
					diagElements[k0]=diagElements[k0]+Ei;
				}
			}
		}
		imp_num=imp_num+1;
	}
	
	std::cout<<"The number of impurities is: "<<imp_num<<std::endl;

	std::ofstream puddle_output("puddle.dat");

	for(int k0=0;k0<diagElements.size();k0++)
	{
		int i0,i1,is,io;
		puddle::IndexToIndexes(	k0, i0, n0 , i1,  n1 ,  io, norb, is, nspin );

		if ( is ==0)
		{
			puddle_output<< i0*A[0][0] + i1*A[1][0] + io*Delta[0]<<" " ;
			puddle_output<< i0*A[0][1] + i1*A[1][1] + io*Delta[1]<<" "  ;
			puddle_output<<diagElements[k0]<<std::endl;
		}
	}
	puddle_output.close();
 return 0;
 }

#endif
