#ifndef KPM_MEMORY_HPP
#define KPM_MEMORY_HPP

#include <cstdlib>
#include <cstring>

#define ALIGN 64

int aligned_malloc( void** p,const size_t required_bytes,const size_t alignment)
{
#ifdef _MKL_H_
	*p= mkl_malloc( required_bytes, alignment );
	if( *p != NULL ) return 0;
	else return -1; 
#else
	return posix_memalign(  p,alignment, required_bytes);
#endif
	
}

void aligned_free(void *p)
{
#ifdef _MKL_H_
    mkl_free(p); p=NULL;
#else
    std::free(p);p=NULL;
#endif
}


void CreateAlignedMemory(const size_t Dim, kpm::complex** memblock )
{
	//RESERVE MEMORY FOR KPM VECTORS
	const size_t
	alignment=64, 
	memsize=( Dim )*sizeof(*memblock[0]);
	
	//Print if one get a status = 0 (ERROR) 
	if ( aligned_malloc( (void**)&memblock[0] ,memsize,alignment)  != 0 )
	{
		std::cerr<<"Unable to allocate enough memory: "<< ((memsize/1024.)/1024.)/1024. <<" GB."<<std::endl;
		std::cerr<<" Aborting simulation"<<std::endl;
		std::exit(-1);
	}
	//Set everything to zero
	memset(memblock[0],0,Dim*sizeof(memblock[0][0])); //Clear storage vector
}


#endif
