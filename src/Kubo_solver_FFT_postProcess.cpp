#include "Kubo_solver_FFT.hpp"
#include "time_station.hpp"

namespace chebyshev{

  
  
Kubo_solver_FFT_postProcess::Kubo_solver_FFT_postProcess(Kubo_solver_FFT& parent_solver): parent_solver_(parent_solver){
  int nump = parent_solver_.nump();
  
  E_points_ = new r_value_t [nump];

  for(int k=0; k<nump;k++)
    E_points_[k] = cos(M_PI * ( 2 * k + 0.5 ) / nump );
  
};



void Kubo_solver_FFT_postProcess::operator()(const value_t final_data[], const value_t r_data[], int r){
  if(parent_solver_.SolverFormula() == KUBO_GREENWOOD)
    Greenwood_postProcess(final_data, r_data, r);
  if(parent_solver_.SolverFormula() == KUBO_BASTIN)
    Bastin_postProcess(final_data, r_data, r);
};

  
void Kubo_solver_FFT_postProcess::eta_CAP_correct(r_value_t E_points[], r_value_t r_data[]){
  int nump = parent_solver_.nump();
  
  for(int e=0; e < nump; e++ )
    r_data[e] *= sin( acos( E_points[e] ) );
}


void Kubo_solver_FFT_postProcess::integration_linqt(const r_value_t E_points[], const r_value_t integrand[], r_value_t data[]){

  int nump = parent_solver_.nump();

  double acc = 0;
  for( int i=0; i < nump-1; i++)
  {
     const double denerg = E_points[i+1]-E_points[i];
     acc +=integrand[i]*denerg;
     data[i] = acc;
  }
}

  


void Kubo_solver_FFT_postProcess::integration(const r_value_t E_points[], const r_value_t integrand[], r_value_t data[]){

  int M     = parent_solver_.num_vectors(),
      nump = parent_solver_.nump();
  
  r_value_t edge = 1.0 - CUTOFF;//parameters_.edge_;

  //#pragma omp parallel for 
  for(int k=0; k<nump - int( M * edge / 4.0 ); k++ ){  //At the very edges of the energy plot the weight function diverges, hence integration should start a little after and a little before the edge.
                                       //The safety factor guarantees that the conductivity is zero way before reaching the edge. In fact, the safety factor should be used to determine
                                      //the number of points to be ignored in the future;
    for(int j=k; j<nump-int(M*edge/4.0); j++ ){//IMPLICIT PARTITION FUNCTION. Energies are decrescent with e (bottom of the band structure is at e=M);
      r_value_t ej  = E_points[j],
	ej1      = E_points[j+1],
	de       = ej-ej1,
        integ    = ( integrand[j+1] + integrand[j] )/2;     
      
      data[k] +=  de * integ;
    }
  }
}


  
void Kubo_solver_FFT_postProcess::partial_integration(const r_value_t E_points[], const r_value_t integrand[], r_value_t data[]){

  int M     = parent_solver_.num_vectors(),
      end   = M * 0.55,
      start = M * 0.45;
  
  //#pragma omp parallel for 
  for(int k=start; k<end; k++ ){  //At the very edges of the energy plot the weight function diverges, hence integration should start a little after and a little before the edge.
                                       //The safety factor guarantees that the conductivity is zero way before reaching the edge. In fact, the safety factor should be used to determine
                                      //the number of points to be ignored in the future;
    for(int j=k; j<end; j++ ){//IMPLICIT PARTITION FUNCTION. Energies are decrescent with e (bottom of the band structure is at e=M);
        r_value_t ej  = E_points[j],
	ej1      = E_points[j+1],
	de       = ej-ej1,
        integ    = ( integrand[j+1] + integrand[j] ) / 2.0;     
      
      data[k] +=  de * integ;
    }
  }
}
     

void Kubo_solver_FFT_postProcess::rearrange_crescent_order( r_value_t* rearranged){//The point choice of the FFT has unconvenient ordering for file saving and integrating; This fixes that.
  int nump = parent_solver_.nump();
  r_value_t original[nump];

  for( int k=0; k < nump; k++ )
    original[k] = rearranged[k];

  
  for( int k=0; k < nump / 2; k++ ){
    rearranged[ 2 * k ]   = original[ k ];
    rearranged[ 2 * k + 1 ] = original[ nump - k - 1 ]; 
  }
  
  for( int k=0; k < nump / 2; k++ ){
    r_value_t tmp = rearranged[ k ]; 
    rearranged[ k ]   = rearranged[ nump-k-1 ];
    rearranged[ nump-k-1 ] = tmp;    
  }  
}


  
void Kubo_solver_FFT_postProcess::Bastin_postProcess(const value_t final_data[], const value_t r_data[],  int r){

  const std::complex<double> im(0,1);  

  
  std::string filename = parent_solver_.OutputFilename();
  
  
  int nump = parent_solver_.nump();
    //R = parameters_.R_,
    //D = parameters_.dis_real_;
  
  int DIM = parent_solver_.Hamiltonian().SystemSize();    

  r_value_t a = parent_solver_.Hamiltonian().HalfWidth(),
    b = parent_solver_.Hamiltonian().BandCenter();
    //    sysSubLength = device_.sysSubLength();
  
  r_value_t omega = DIM/( a * a );//* sysSubLength * sysSubLength );//Dimensional and normalizing constant
  
  //r_value_t tmp, max=0, av=0;

  r_value_t
    rearranged_E_points[nump],
    integrand[nump],
    rvec_integrand[nump],
    partial_result[nump],
    rvec_partial_result[nump];

  for(int e=0; e<nump; e++){
    rearranged_E_points[e] = E_points_[e];
    integrand[e]           = 0.0;
    rvec_integrand[e]      = 0.0;
    partial_result[e]      = 0.0;
    rvec_partial_result[e] = 0.0;
  }   

  rearrange_crescent_order(rearranged_E_points);
  /*When introducing a const. eta with modified polynomials, the result is equals to that of a
  simulation with regular polynomials and an variable eta_{var}=eta*sin(acos(E)). The following
  heuristical correction greatly improves the result far from the CNP to match that of the
  desired regular polys and const. eta.*/
  


  /*//R convergence analysis  
  for(int e = 0; e < 2 * nump; e++){

    tmp = real( final_data[ e ] );
    final_data_[e] = ( final_data_ [ e ] * ( r - 1.0 ) +  r_data_[ e ] ) / r;

    if( r > 1 ){

      tmp = std::abs( ( final_data_ [ e ] - tmp ) / tmp) ;
      if(tmp > max)
        max = tmp;

      av += tmp / nump ;
    }
  }

  
  if(r>1){
    conv_R_[ 2 * ( r - 1 ) ]   = max;
    conv_R_[ 2 * ( r - 1 ) + 1 ] = av;
  }
  */

    

  //Keeping just the real part of E*p(E)+im*sqrt(1-E^2)*w(E) yields the Kubo-Bastin integrand:
  for(int k = 0; k < nump; k++){
    integrand[k]  = E_points_[k] * real( final_data[ k ] ) - ( sqrt(1.0 - E_points_[ k ] * E_points_[ k ] ) * imag( final_data[ k + nump ] ) );
    integrand[k] *= 1.0 / pow( (1.0 - E_points_[k]  * E_points_[k] ), 2.0);
    integrand[k] *=  omega / ( M_PI ); 

    rvec_integrand[k]  = E_points_[k] * real( r_data[ k ] ) - ( sqrt(1.0 - E_points_[ k ] * E_points_[ k ] ) * real( r_data[ k + nump ] ) );
    rvec_integrand[k] *= 1.0 / pow( (1.0 - E_points_[k]  * E_points_[k] ), 2.0);
    rvec_integrand[k] *=  omega / ( M_PI ); 
  }

  rearrange_crescent_order(integrand);
  rearrange_crescent_order(rvec_integrand);

  


  time_station time_integration;
      
  integration_linqt(rearranged_E_points, rvec_integrand, rvec_partial_result);  
  integration_linqt(rearranged_E_points, integrand, partial_result);    

  time_integration.stop("       Integration time:           ");

  

  /*
  if( parameters_.eta_!=0 ){
    eta_CAP_correct(rearranged_E_points, partial_result);
    eta_CAP_correct(rearranged_E_points, rvec_partial_result);
  }
  */



  
  /*
  std::ofstream dataR;
  dataR.open(run_dir+"vecs/r"+std::to_string(r)+"_"+filename);

  for(int e=0;e<nump;e++)  
    dataR<< a * rearranged_E_points[e] - b<<"  "<< rvec_partial_result [e] <<std::endl;

  dataR.close();
  */


  
  std::ofstream dataP;
  dataP.open(filename);

  for(int e=0;e<nump;e++)  
    dataP<< a * rearranged_E_points[e] - b<<"  "<<  partial_result [e] <<std::endl;

  dataP.close();



  
  /*
  std::ofstream data;
  data.open(run_dir+"conv_R_"+filename);

  for(int r = 1; r < D * R; r++)  
    data<< r <<"  "<< conv_R_[ 2 * ( r - 1 ) ]<<"  "<< conv_R_[ 2 * ( r - 1 ) + 1 ] <<std::endl;

  data.close();
  */

  
  
  
  std::ofstream data2;
  data2.open(filename+"_integrand");

  for(int e=0;e<nump;e++)  
    data2<< a * rearranged_E_points[e] - b<<"  "<<  integrand[e] <<std::endl;
  
  data2.close();


  //  plot_data(run_dir,filename);

}





void Kubo_solver_FFT_postProcess::Greenwood_postProcess(const value_t final_data[], const value_t r_data[], int r){

  const std::complex<double> im(0,1);  

  std::string filename = parent_solver_.OutputFilename();

  
  int nump = parent_solver_.nump();
  //    R = parameters_.R_,
  //  D = parameters_.dis_real_;
  
  int DIM = parent_solver_.Hamiltonian().SystemSize();    

  r_value_t a = parent_solver_.Hamiltonian().HalfWidth(),
    b = parent_solver_.Hamiltonian().BandCenter();
  //    sysSubLength = device_.sysSubLength();
  
  r_value_t omega = DIM/( a * a );//* sysSubLength * sysSubLength );//Dimensional and normalizing constant
  
  //  r_value_t tmp, max=0, av=0;

  r_value_t
    rearranged_E_points[nump],
    partial_result[nump],
    rvec_partial_result[nump];

  
  for(int e=0; e<nump; e++){
    rearranged_E_points[e] =  E_points_[e] ;
    partial_result[e] = 0.0;
    rvec_partial_result[e] = 0.0;
  }
  rearrange_crescent_order(rearranged_E_points);
  /*When introducing a const. eta with modified polynomials, the result is equals to that of a
  simulation with regular polynomials and an variable eta_{var}=eta*sin(acos(E)). The following
  heuristical correction greatly improves the result far from the CNP to match that of the
  desired regular polys and const. eta.*/


  /*
  for(int e = 0; e < nump; e++){

    tmp = real( final_data_[ e ] );
    final_data_[e] = ( final_data_ [ e ] * ( r - 1.0) + r_data_[ e ] ) / r;

    if( r > 1 ){

      tmp = std::abs( ( final_data_ [ e ] - tmp ) / tmp) ;
      if(tmp > max)
        max = tmp;

      av += tmp / nump ;
    }
  }
  
  
  if(r>1){
    conv_R_[ 2 * (r-1) ]   = max;
    conv_R_[ 2 * (r-1)+1 ] = av;
  }
  */

  
  for(int k=0; k < nump; k++){
    rvec_partial_result[k] = 2.0 * omega * real( r_data[k] )     / (1.0 - E_points_[k] * E_points_[k] ) / ( 2 * M_PI );
    partial_result[k]      = 2.0 * omega * real( final_data[k] ) / (1.0 - E_points_[k] * E_points_[k] ) / ( 2 * M_PI );
  }

  rearrange_crescent_order(partial_result);
  rearrange_crescent_order(rvec_partial_result);


  /*
  if( parameters_.eta_ != 0 ){
    eta_CAP_correct(rearranged_E_points, partial_result);
    eta_CAP_correct(rearranged_E_points, rvec_partial_result);
  }
  */


  
  /*
  std::ofstream dataR;
  dataR.open(run_dir+"vecs/r"+std::to_string(r)+"_"+filename);

  for(int e=0;e<nump;e++)  
    dataR<< a * rearranged_E_points[e] - b<<"  "<<  rvec_partial_result [e] <<std::endl;

  dataR.close();
  */


  
  std::ofstream dataP;
  dataP.open(filename);

  for(int e=0;e<nump;e++)  
    dataP<<  a * rearranged_E_points[e]-b<<"  "<< partial_result [e] <<std::endl;

  dataP.close();



  
  /*
  std::ofstream data;
  data.open(run_dir+"conv_R_"+filename);

  for(int r = 1; r < D * R; r++)  
    data<< r <<"  "<< conv_R_[ 2 * ( r - 1 ) ]<<"  "<< conv_R_[ 2 * ( r - 1 ) + 1 ] <<std::endl;

  data.close();
  */

  

  //  plot_data(run_dir,filename);  

}








void Kubo_solver_FFT_postProcess::plot_data(const std::string& run_dir, const std::string& filename){
        //VIEW commands
  
     std::string exestring=
         "gnuplot<<EOF                                               \n"
         "set encoding utf8                                          \n"
         "set terminal pngcairo enhanced                             \n"

         "unset key  \n"

         "set output '"+run_dir+filename+".png'                \n"

         "set xlabel 'E[eV]'                                               \n"
         "set ylabel  'G [2e^2/h]'                                           \n"
         
        "plot '"+run_dir+"currentResult_"+filename+"' using 1:2 w p ls 7 ps 0.25 lc 2;  \n"
         "EOF";
     
      char exeChar[exestring.size() + 1];
      strcpy(exeChar, exestring.c_str());    
      if(system(exeChar)){};


}

};
