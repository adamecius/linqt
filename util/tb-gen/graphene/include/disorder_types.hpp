#ifndef RADIAL_DISORDER
#define RADIAL_DISORDER

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

int radial_disorder(	
		const int n0, const int n1, 
		const int norb, const int nspin,
		Real (*maxHeight)(int, int,int,int),
		Real concen,
		Real range,	
		Real (*PotProfile)(int, int,int,int,Real,Real)
		std::vector<Real>& diagElements, Real seed)
{
	srand((int)seed);
	const int dim=n0*n1*norb*nspin;
	const double epsilon=std::numeric_limits<double>::epsilon();


	//Check if the height or the concentration are not zero
	if(maxHeight <= epsilon || concen  <= epsilon)
		return 0;

	//Estimate the cutoff indexes		
	const Real 
	L0=sqrt( pow(A[0][0],2.) + pow(A[0][1],2.) ) ,
	L1=sqrt( pow(A[1][0],2.) + pow(A[1][1],2.) ) ;

	//Evaluate how many unit cells would be necessary 
	//using range. In general range is not a
	// hard boundary, therefore we include a range tolerance
	//parameter
	const Real range_tol=20; 
	const Real CellCutOff = range_tolrange_tol*Range ; 

	//Dividing by the max distance gives you the 
	// the number of cells
	const int 
	i0c=CellCutOff/L0  ,  
	i1c=CellCutOff/L1  ;  
	std::cout<<"The disorder will be evaluated using: "<<std::endl;
	std::cout<<"super cells : "<<i0c<<" "<<i1c<<std::endl;

	//Incluiding disorder
	//We go through all the possible positions
	int imp_num=0; //used to count the number of impurities
	for(int ip0=0; ip0 < n0 ; ip0++ )
	for(int ip1=0; ip1 < n1  ; ip1++ )
	for(int ipo=0; ipo < norb ; ipo++ )
	if(  (double)rand()/(double)RAND_MAX < concen || (ipo0== n0/2&& ip1==n1/2 && ipo==0) ) //Select a center
	{
		const double 
		Umax = maxHeight(ip0,ip1,ipo,0);
		
		//The we include the effect of disorder in the 
		//neighboring cells
		for(int i0=-i0c; i0 <= i0c ; i0++ )
		for(int i1=-i1c; i1 <= i1c ; i1++ )
		for(int IO=0; IO < norb ; IO++ )
		{
			const int 
			I0 = ip0+i0,
			I1 = ip1+i1;

			//Compute the distance from the center
			const double
			rdiff[2] ={ (ip0-I0)*A[0][0] + (ip1-I1)*A[1][0]+ (IO-ipo)*Delta[0],
						(ip0-I0)*A[0][1] + (ip1-I1)*A[1][1]+ (IO-ipo)*Delta[1]
						},
			dist = 	sqrt( rdiff[0]*rdiff[0] + rdiff[1]*rdiff[1] );
			
			//If the distance is smaller than CellCutOff we include the disordr
			if( dist < CellCutOff*0.5 )
			{
				//because is scalar, its the same for both spins
				for( int IS=0; IS<nspin; IS++)
				{
					const double ss= 1.0 - 2.0*IS ;
					const double
					Ei= PotProfile( I0,I1,IO,IS,Umax,dist);
					const int
					k0= puddle::IndexesToIndex(	I0, n0,
												I1, n1,
												IO, norb,
												IS, nspin );
					diagElements[k0]+=Ei;
				}
			}
		}
		imp_num=imp_num+1;
	}
 return 0;
 }



class PuddleDisorder( 	
{
	public:
	
	PuddleDisorder( Real _height, Real _concen, Real _range):
		height(_height), concen(_concen), range(_range)
		{};
		
	Real maxHeight( const int i0, const int i1, 
					const int io, const int is)
	{
		return height*(2.*( (float)rand()/(float)RAND_MAX ) - 1.);
	}

	Real PotProfile(const int i0, const int i1, 
					const int io, const int is,
					Real Umax,Real dist)
	{
		return Umax*exp( -0.5*dist*dist/range/range );
	}

	Real height, concen, range;
};
						
class PointDisorder
{
	public:
	
	PointDisorder( Real _height, Real _concen, Real _range):
		height(_height), concen(_concen), range(_range)
		{};
		
	Real maxHeight( const int i0, const int i1, 
					const int io, const int is)
	{
		return height*(2.*( (float)rand()/(float)RAND_MAX ) - 1.);
	}

	Real PotProfile(const int i0, const int i1, 
					const int io, const int is,
					Real Umax,Real dist)
	{
		if( dist == 0.0)
			return Umax;
		else
			return 0.0;
	}

	Real height, concen, range;
};

class StaggeredDiskDisorder
{
	public:
	
	StaggeredDiskDisorder( Real _height, Real _concen, Real _range):
		height(_height), concen(_concen), range(_range)
		{};
		
	Real maxHeight( const int i0, const int i1, 
					const int io, const int is)
	{
		return height*( 1. - 2.*io .);
	}

	Real PotProfile(const int i0, const int i1, 
					const int io, const int is,
					Real Umax,Real dist)
	{
		if( dis< range 0.0)
			return Umax;
		else
			return 0.0;
	}

	Real height, concen, range;
};

#endif
