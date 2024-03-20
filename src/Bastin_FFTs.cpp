#include "fftw_wrapper.hpp"
#include "Kubo_solver_FFT.hpp"

namespace chebyshev{

void Kubo_solver_FFT::Bastin_FFTs(Vectors_sliced &chebVecL, Vectors_sliced &chebVecR, value_t r_data[], int s){

  const std::complex<double> im(0,1);

  
  size_t size = section_size_;
  
  if( s != num_sections_-1 )
    size -= Hamiltonian_dummyMoms_.SystemSize() % num_sections_;


  
  value_t pre_factors [ nump_ ];      
      
  for(int m = 0; m < M_; m++)
    pre_factors[m]  = value_t( ( 2 - ( m == 0 ) ) * Jackson_kernel(m, M_) ) * std::polar( 1.0, M_PI * m / ( 2.0 * nump_ ) ) ;




#pragma omp parallel 
  {
    int id,  Nthrds, l_start, l_end;
    id      = omp_get_thread_num();
    Nthrds  = omp_get_num_threads();
    l_start = id * size / Nthrds;
    l_end   = (id+1) * size / Nthrds;

    if (id == Nthrds-1)
      l_end = size;
  
  
  //8 plans + 14 [M] sized vectors per thread. Can it be reduced to 2 plans and some 4 vectors?

    out_of_place_dft
      re_bras( nump_, BACKWARD ),
      im_bras( nump_, BACKWARD ),
      
      re_kets( nump_, BACKWARD ),
      im_kets( nump_, BACKWARD ),

     
      re_D_bras( nump_, BACKWARD ),
      im_D_bras( nump_, BACKWARD ),
      
      re_D_kets( nump_, BACKWARD ),
      im_D_kets( nump_, BACKWARD );


   std::complex<double> bra, ket, //variables to access values within bras/kets         
                        p[nump_], w[nump_]; //Dot product partial results;

   for(int k = 0; k < nump_; k++){
     p[k] = 0;
     w[k] = 0;
   }
    

    //FFTW plans. All plans are in-place backward 1D FFTS of the entry variables; These functions are no thread safe for some reason, hence the # pragma critical
# pragma omp critical
   {
     re_bras.create();
     im_bras.create();
     re_kets.create();
     im_kets.create();
     
     re_D_bras.create();
     im_D_bras.create();
     re_D_kets.create();
     im_D_kets.create();
   }

   for(int l = l_start; l < l_end; l++){
     for(int m = 0; m < M_; m++){
	//The following casting allows to simplify all template instantiations: converts all to double.
	//All operations here are performed in double precision.
	//I also don`t think fftw supports long double as of right now.
	//The row-wise accessed of the col-major matrices are some of the longest operations in this loop!!

       std::complex<r_value_t> im_m = static_cast<std::complex<r_value_t>> (m);
	
       bra = chebVecL.Vector(m)[l];

       //Bra Greens functions FFT inputs:
       re_bras.input()[m] = pre_factors[m] * real( bra );
       im_bras.input()[m] = pre_factors[m] * imag( bra );
	  
	
       //The same inputs, multiplied by m, serve as inputs for the derivative FFT: 
       re_D_bras.input()[m] = im_m * re_bras.input()[m];
       im_D_bras.input()[m] = im_m * im_bras.input()[m];
	
       

	
       ket = chebVecR.Vector(m)[l];

       //Bra Greens functions FFT inputs:
       re_kets.input()[m] = pre_factors[m] * real( ket );
       im_kets.input()[m] = pre_factors[m] * imag( ket );
	  
	
       //The same inputs, multiplied by m, serve as inputs for the derivative FFT: 
       re_D_kets.input()[m] = im_m * re_kets.input()[ m ];
       im_D_kets.input()[m] = im_m * im_kets.input()[ m ];	
     }

     re_bras.execute();
     im_bras.execute();
     re_kets.execute();
     im_kets.execute();  

     re_D_bras.execute();
     im_D_bras.execute();
     re_D_kets.execute();
     im_D_kets.execute();
  



     //As mentioned previously, variable change as k=2j; This loop updates the partial results of the dot products p(E) and w(E);
     //This also performs the conjugation from bra.cdot(ket) of the complex random vector dot product.
     for(int j = 0; j < nump_/2; j++){
       //Here: p(k) += Re(G(k)) * G(k) + G(k) * Re(G(k)).
       p[j] += ( real( re_bras(j) ) - im * real( im_bras(j) ) ) *  //Re(G(k)) *
	         (       re_kets(j)   + im *       im_kets(j) ) + //G(k)+
 	  
	         ( conj( re_bras(j) ) - im * conj( im_bras(j) ) ) * //G(k)
	         ( real( re_kets(j) ) + im * real( im_kets(j) ) );   //Re(G(k))
	

       //Here: w(k) += (dG(k)) * Re(G(k)) - Re(G(k)) * (dG(k))
       w[j] += ( conj( re_D_bras(j) ) - im * conj( im_D_bras(j) ) )  * //dG(k)
		 ( real( re_kets(j)   ) + im * real( im_kets(j) ) ) - //Re(G(k))
	  
		 ( real( re_bras(j) )   - im * real( im_bras(j) ) )  *  //Re(G(k))
	         ( re_D_kets(j)         + im * im_D_kets(j) ); //dG(k)
     }

     
     for(int j = nump_ / 2; j < nump_; j++){
       //Here: p(k) += Re(G(k)) * G(k) + G(k) * Re(G(k)).
       p[j] += ( real( re_bras(j) ) - im * real( im_bras(j) ) ) *  //Re(G(k)) *

	       ( conj(re_kets(j) )   + im *       conj( im_kets(j) ) ) + //G(k)+
 	  
	         ( re_bras(j)   - im * im_bras(j)  ) * //G(k)
	         ( real( re_kets(j) ) + im * real( im_kets(j) ) );   //Re(G(k))
	

         //Here: w(k) += (dG(k)) * Re(G(k)) - Re(G(k)) * (dG(k))
         w[j] += (  re_D_bras(j)  - im *  im_D_bras(j)  )  * //dG(k)
		 ( real( re_kets(j)   ) + im * real( im_kets(j) ) ) - //Re(G(k))
	  
		 ( real( re_bras(j) )   - im * real( im_bras(j) ) )  *  //Re(G(k))
	         ( conj( re_D_kets(j) )         + im * conj ( im_D_kets(j) ) ); //dG(k)

     }
   }      
    
# pragma omp critical
   {
     for(int k = 0; k < nump_; k++){
       r_data[ k ]         += p[ k ];
       r_data[ k + nump_ ] += w[ k ];	
     }
   } 
 }
}
  
};
