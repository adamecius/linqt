#ifndef ANDERSON_DISORDER
#define ANDERSON_DISORDER

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


namespace Anderson
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

int AndersonDisorder(	
			const int n0, const int n1, 
			const int norb, const int nspin,
			const Real Wmax,
			std::vector<Real>& diagElements, Real seed)
{
	srand((int)seed);
	const int dim=n0*n1*norb*nspin;
	const double epsilon=std::numeric_limits<double>::epsilon();

	// first we go through all the possible positions
	for(int i0=0; i0 < n0  ; i0++ )
	for(int i1=0; i1 < n1  ; i1++ )
	for(int io=0; io < norb; io++ )
	{
		const Real W= Wmax*(( (Real)rand()/(Real)RAND_MAX ) - 0.5); ;
		for( int is=0; is<nspin;is++)
		{
			int k0= Anderson::IndexesToIndex( i0, n0, i1, n1, io, norb, is, nspin);
			diagElements[k0]=diagElements[k0]+W;
		}
	}
 return 0;
 }

#endif
