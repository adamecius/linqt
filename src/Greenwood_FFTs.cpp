#include "fftw_wrapper.hpp"
#include "Kubo_solver_FFT.hpp"

namespace chebyshev{

void Kubo_solver_FFT::Greenwood_FFTs(Vectors_sliced &chebVecL, Vectors_sliced &chebVecR, value_t r_data[], int s){

  size_t size = section_size_;
  
  if( s != num_sections_-1 )
    size -= Hamiltonian_dummyMoms_.SystemSize() % num_sections_;


 
  value_t pre_factors [ M_ ];

  for(int m = 0; m < M_; m++)
    pre_factors[m]  = value_t( ( 2 - ( m == 0 ) ) * Jackson_kernel(m, M_) ) * std::polar( 1.0, M_PI * m / ( 2.0 * nump_ ) ) ;


  
  const std::complex<double> im(0,1);

  
#pragma omp parallel 
  {
    int id,  Nthrds, l_start, l_end;
    id      = omp_get_thread_num();
    Nthrds  = omp_get_num_threads();
    l_start = id * size / Nthrds;
    l_end   = (id+1) * size / Nthrds;
    
    if (id == Nthrds-1) l_end = size;


    
    value_t thread_data [nump_];

    for(int k=0;k<nump_;k++)
      thread_data[k]=0;

    
    out_of_place_dft
      re_bras( nump_, BACKWARD ),
      im_bras( nump_, BACKWARD ),
      
      re_kets( nump_, BACKWARD ),
      im_kets( nump_, BACKWARD );
     
    	
# pragma omp critical
    {
      re_bras.create();
      im_bras.create();
      re_kets.create();
      im_kets.create();
    }

    for(int l = l_start; l < l_end; l++){

      for( int m = 0; m < M_; m++ ){
	re_bras.input()[m] = pre_factors[m] * real( chebVecL.Vector(m)[l] );
        im_bras.input()[m] = pre_factors[m] * imag( chebVecL.Vector(m)[l] );

	re_kets.input()[m] = pre_factors[m] * real( chebVecR.Vector(m)[l] );
        im_kets.input()[m] = pre_factors[m] * imag( chebVecR.Vector(m)[l] ); 
      }

      
      re_bras.execute();
      im_bras.execute();
      re_kets.execute();
      im_kets.execute();


      for(int k=0; k<nump_; k++ ){
        thread_data[k] += real (
				 ( real( re_bras(k) ) - im * real( im_bras(k) ) ) *
				 ( real( re_kets(k) ) + im * real( im_kets(k) ) )
			       );

      }
    }

    # pragma omp critical
    {
      for(int k = 0; k < nump_; k++)
	r_data[k] += thread_data[k];
    }
  }
}


}
