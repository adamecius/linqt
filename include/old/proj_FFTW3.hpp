#ifndef PROJ_FFTW3_HPP
#define PROJ_FFTW3_HPP


#include "types_definitions.hpp"
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>
#include <fftw3.h>

class BZProjection
{
	
	public: 
	//Constructor
	BZProjection( ):rank(3)
	{
		in_= NULL;
		out_= NULL;
		filterIdx_=NULL;
	};

	
	
	//Methods
	void CleanUp()
	{
		fftw_destroy_plan(planForw_);
		fftw_destroy_plan(planBack_);
		
		if(in_!= NULL) fftw_free(in_); 
		if(out_!= NULL) fftw_free(out_);
		if(filterIdx_!=NULL) free(filterIdx_);
		fftw_cleanup();
	}

	void StartProjection( std::string _projLabel)
					{
						projLabel_=_projLabel;
						std::string projFileName_="operators/"+projLabel_+".KPROJ";
						std::ifstream projFile;
						projFile.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
			
						//Try to open the file projFileName_  
						
						try  {	projFile.open ( projFileName_.c_str() );	}
						catch (std::ifstream::failure e)
						{
							projFile.close();
							std::cerr 	<<"Exception:"<<e.what()
										<<"\n while opening config file: "
										<< projFileName_<<"."<<std::endl;
						}
						
						int orbNum, maxSpin, orbPerCell, numCells;
						projFile>>n[0]>>n[1]>>n[2]>> orbNum>> maxSpin;	
						numCells 	= n[0]*n[1]*n[2];
						orbPerCell 	= orbNum*maxSpin;
						dim_		= numCells*orbPerCell;
						
						//Initialize the input and output state vectors
						in_  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Dim() );
						out_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Dim() );
						filterIdx_ = (double*) malloc( sizeof(double) * Dim() ); 

						//Define the dimension for each FFT
						const int howmaby = orbPerCell; //We need to perform a FFT for each orbital
						int* nembed = n;	//Size of the sub-array
						int dist = 1;//distance between arrays
						int stride= orbPerCell;	//The array is not continous 
						std::string planLabel=(projLabel_+".FFTF");
					
						//Choose the best plans for the particular system

						bool foundPlan = fftw_import_wisdom_from_filename(planLabel.c_str()) ;

						planForw_ 	=fftw_plan_many_dft(rank,n, howmaby, 
													 in_ , nembed, stride, dist,
													 out_, nembed, stride, dist,
													 FFTW_FORWARD, FFTW_PATIENT);

						planBack_ 	=fftw_plan_many_dft(rank,n, howmaby, 
													 out_ , nembed, stride, dist,
													 in_  , nembed, stride, dist,
													 FFTW_BACKWARD, FFTW_PATIENT);
						if(!foundPlan)
							fftw_export_wisdom_to_filename(planLabel.c_str());
					

						projFile.exceptions ( std::ifstream::badbit );
						while ( !projFile.eof() )
						{
							int n0; 
							double Idx;
							projFile>>n0>>Idx;	
							for( int orb= 0 ; orb< orbPerCell&& !projFile.eof() ; orb ++ )
								filterIdx_[n0*orbPerCell+orb]=Idx; 
						}
						projFile.close();						

						//for testing purposes
//						for( int n0=0; n0< numCells; n0++)
//						for( int orb= 0 ; orb< orbPerCell ; orb ++ )
//							filterIdx_[n0*orbPerCell+orb]=1.0; 
												
				};

	
	void DirectFourierTransform( int size, std::complex<double>*   X	)
	{
		memcpy( &in_[0] , &X[0]   , Dim()*sizeof( fftw_complex ) );//Pass the vector to init
		fftw_execute( planForw_ ); //Transfor the vector from real to momentum space
		memcpy( &X[0]   , &out_[0], Dim()*sizeof( fftw_complex ) );//Pass Final result to X
	}

	void InverseFourierTransform( int size,std::complex<double>*   X)
	{
		memcpy( &out_[0] , &X[0]   , Dim()*sizeof( fftw_complex ) );//Pass the vector to init
		fftw_execute( planBack_ ); //Transfor the vector from real to momentum space
		memcpy( &X[0]   , &in_[0], Dim()*sizeof( fftw_complex ) );//Pass Final result to X
	}
			
	//Returns the total dimension of the system
	const int Dim()
	{
		return dim_;
	}

	void ProjectVector( int size,std::complex<double>*   X)
	{
		if (size != Dim() )
			std::cerr<<"The dimension of the input vetor in "
					 <<"ProjectVector method of Projection class does "
					 <<"not mathc internal dimension"<<std::endl;
	
		memcpy( &in_[0], &X[0]   , Dim()*sizeof( fftw_complex ) );//Pass the vector to init		
		fftw_execute( planForw_ ); //Transfor the vector from real to momentum space
		for( int i=0;i < Dim(); i++) //Apply projection
		{	
			(out_[i])[0] = (out_[i])[0]*filterIdx_[i];
			(out_[i])[1] = (out_[i])[1]*filterIdx_[i];
		}		
		fftw_execute( planBack_ ); //Transfor the vector from momentum to real space
		memcpy( &X[0], &in_[0] , Dim()*sizeof( fftw_complex ) );
 
		//Pass Final result to X and normalize
		double numCells = n[0]*n[1]*n[2];
		for(int i=0; i< Dim() ;i++)
			X[i]/=numCells;


	}
	
	private: 
	int rank ;
	int dim_, n[3];
	fftw_complex *in_,*out_;
	double* filterIdx_;
	fftw_plan planForw_, planBack_ ;
	std::string projLabel_;

};
#endif
