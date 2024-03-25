#include "Kubo_solver_FFT.hpp"
#include "time_station.hpp"



namespace chebyshev 
{



  
void Kubo_solver_FFT::compute( SparseMatrixType &OPL, SparseMatrixType &OPR,  qstates::generator& gen ){

  time_station solver_station;


  const size_t DIM = chebVecL_.SystemSize();
  Kubo_solver_FFT_postProcess postProcess(*this);
  section_size_ = DIM / num_sections_ +
                  DIM % num_sections_; //We are assuming */* is bigger than *%*. If not, i'd increase/decrease num_parts until it is.
  


  
  
  time_station allocation_time;

  //Allocate the memory
  
  value_t
    r_data[2*nump_],
    final_data[2*nump_];	

  allocate(chebVecL_, chebVecR_);

  
  allocation_time.stop("\n \nAllocation time:            ");





  
  
  gen.SystemSize(DIM);

  while( gen.getQuantumState() )  {
    time_station randVec_time;
    std::cout<<"Computing with ID: "<<gen.count<<" states" <<std::endl;

    //SELECT RUNNING TYPE
    polynomial_recursion(gen.State(),gen.State(), OPL, OPR, chebVecL_,chebVecR_, r_data);

    randVec_time.stop("       Total RandVec time:         ");
    std::cout<<std::endl;



     
     time_station time_postProcess;
     
     update_data(final_data, r_data, gen.count);
     postProcess(final_data, r_data, gen.count);
     reset_data(r_data);
     
     time_postProcess.stop( "       Post-processing time:       ");


  }
	
  solver_station.stop("Total case execution time:              ");
}


  

  

void Kubo_solver_FFT::polynomial_recursion(const vector_t& PhiR, const vector_t& PhiL,
				SparseMatrixType &OPL,
				SparseMatrixType &OPR,  
				Vectors_sliced &chebVecL,
				Vectors_sliced &chebVecR,
				value_t r_data[]){


        int total_time_csrmv = 0,
            total_time_FFTs  = 0;
	



        for(int s=0; s < num_sections_; s++){

	  std::cout<< "   -Section: "<<s+1<<"/"<<num_sections_<<std::endl;


	  
	  time_station csrmv_time_kets;
	  
	  chebVecL.SetInitVectors( OPL, PhiL );
	  chebVecL.IterateAllSliced(s);

	  csrmv_time_kets.stop_add( &total_time_csrmv, "           Kets cycle time:            ");


	  

	  time_station csrmv_time_bras;  

	  chebVecR.SetInitVectors(  PhiR );
	  chebVecR.IterateAllSliced(s);
	  chebVecR.MultiplySliced( OPR, s );
		
	  csrmv_time_bras.stop_add( &total_time_csrmv,  "           Bras cycle time:            ");

	
  
  
	  time_station FFTs_time;

	  if( sym_formula_ == KUBO_GREENWOOD )
	    Greenwood_FFTs(chebVecL, chebVecR, r_data, s);

          if( sym_formula_ == KUBO_BASTIN )
	    Bastin_FFTs   (chebVecL, chebVecR, r_data, s);	

	  FFTs_time.stop_add( &total_time_FFTs, "           FFT operations time:        ");

	
	}

      
        time_station total_CSRMV(total_time_csrmv, "\n       Total CSRMV time:           ");
        time_station total_FFTs(total_time_FFTs, "       Total FFTs time:            ");
  }  






  
void Kubo_solver_FFT::allocate(Vectors_sliced& chebVecL, Vectors_sliced& chebVecR){
  const size_t DIM = Hamiltonian_dummyMoms_.SystemSize();

  //This operation is memory intensive
  std::cout<<"Initializing chevVecL"<<std::endl;
  chebVecL.CreateVectorSet( );
  std::cout<<"Initialize chevVecR"<<std::endl;
  chebVecR.CreateVectorSet( );


  
  double buffer_mem    = double( 2 * M_ * section_size_ * sizeof(value_t) ) / double( 1000000000 ),
         recursion_mem = double( ( 3 * DIM  ) * sizeof(value_t) )/ double( 1000000000 ),
         FFT_mem       = 0.0,
         OP_mem = 0, //To be defined;
         Total = 0.0;

  
  if(sym_formula_ == KUBO_GREENWOOD)
    FFT_mem = double( ( 1 + /*num_threads() */ ( 8 + 1 ) ) * nump_ * sizeof(value_t) ) / double( 1000000000 );
  if(sym_formula_ == KUBO_BASTIN)
    FFT_mem = double( ( 1 + /*omp_get_num_threads() */ ( 16 + 1 ) ) * nump_ * sizeof(value_t) ) / double( 1000000000 );

  Total = buffer_mem + OP_mem + recursion_mem + FFT_mem;

  
  std::cout<<std::endl;
  std::cout<<"Expected memory cost breakdown:"<<std::endl;
  std::cout<<"   Chebyshev buffers:    "<< buffer_mem<<" GBs"<<std::endl;  
  std::cout<<"   Operators memory:     "<< OP_mem<<" GBs"<<std::endl;  
  std::cout<<"   Recursion vectors:    "<<  recursion_mem <<" GBs"<<std::endl;
  std::cout<<"   FFT auxiliary lines:  "<<  FFT_mem <<" GBs"<<std::endl<<std::endl;   
  std::cout<<"TOTAL:  "<<  Total<<" GBs"<<std::endl<<std::endl;

}




  //Post-processing




  
};
