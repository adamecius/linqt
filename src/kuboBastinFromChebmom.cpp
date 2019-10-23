
// Used for OPENMP functions
#include <omp.h>

// Parsing library
#include <libconfig.h++>

// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/

// MKL LIBRARIES
#define MKL_Complex16 std::complex<double>
#include "mathimf.h"
#include "mkl.h"
#include "mkl_spblas.h"
typedef std::complex<double> complex;
typedef MKL_INT integer;


//Discrete Cosine Transform
#include <fftw3.h>

//CUSTOM LIBRARIES
#include "kernel_functions.hpp"
#include "algebra_functions.hpp"
#include "bessel_int_coeff.hpp"
#include "time_handler.hpp"

const int INTEGER_NOT_FOUND = 0 ;
const double kb    = 0.086173324; 	// boltlzmann constant meV/K
const double hbar  = 0.6582119624; 	// reduced  planck constant eV.fs


double fermi_dirac( const double x )
{
	if( x > 40. )
		return 0.0;
	
	return 1.0/( exp(x) + 1.0 );  	
};

template <typename T> 
std::string to_string( T x ) 
{
	std::stringstream ss;
	ss << x;
	return ss.str();
}


struct chebMom
{
	chebMom():numMom1(0),numMom0(0){};

	chebMom(const int m0,const int m1):numMom1(m1),numMom0(m0),mu( std::vector<complex>(numMom1*numMom0, 0.0) ){};
	 
	complex& operator()(const int m0,const int m1)
	{
		return mu[ m0*numMom1 + m1 ];
	}
	int numMom1,numMom0;
	std::vector<complex> mu;
};


int main(int argc, char *argv[])
{	
	libconfig::Config cfg;	//Class use to handle files
	cfg.readFile(argv[1]); 	// Read the file. If there is an error, report it and exit.
	
	const
	std::string 
	LABEL = cfg.lookup("SystemName"),
	S_OPR = cfg.lookup("Operator");

	int maxNumMom,maxNumTMom,numRV;
	double xbound = 0.9; //The typical normalization is from (-0.9, 0.9)
	cfg.lookupValue("xbound" ,xbound); std::cout<<"The bound is"<<xbound<<std::endl;
	cfg.lookupValue("NumberOfMoments", maxNumMom );
	cfg.lookupValue("NumberOfTMoments", maxNumTMom );
	cfg.lookupValue("NumberOfStateVec", numRV );

	//Use the data read from input file to determine the prefix of all files
	std::string prefix  ="NonEqOp"+S_OPR+LABEL+"KPM_M"+to_string(maxNumMom)+"x"+to_string(maxNumTMom)+"RV"+to_string(numRV);
	std::string momfilename =prefix+".chebmom2D" ;

	//Open file, read first the number of moments by reading the
	//index of the last element
	std::cout<<std::endl<<"Reading moments from "<<momfilename<<std::endl;
	std::ifstream momfile( momfilename.c_str() );
	//Read the halfwidth and the dimension at the end of the
	//moment file
	double HalfWidth,BandCenter; int dim; 
	momfile>>dim>>HalfWidth>>BandCenter; HalfWidth=HalfWidth/2.0;//in the file what you have is the bandwidth
	std::cout<<"Using (HalfWidth,BandCenter)="<<HalfWidth<<" "<<BandCenter<<std::endl;
	//Read the dimension of the moments and 
	int numMoms0,numMoms1;
	momfile>>numMoms0>>numMoms1;
	//Check if the number of moments in the file is the same as
	//what expected from the configu file
	if( numMoms0 != maxNumMom || numMoms1 != maxNumTMom )
	{
		std::cerr	<<std::endl<<"WARNING!"<<std::endl
					<<" The maximum number of moments obtained from the momfile: "
					<<numMoms0<<","<<numMoms1<<std::endl
					<<"do not correspond with the maximum number of moments from the config file: "<<maxNumMom<<","<<maxNumTMom<<std::endl;
		return -1;
	}
	//If they are, read all the moments
 	chebMom mu(numMoms0,numMoms1);
	for( int m0 = 0 ; m0 < numMoms0 ; m0++)	//energy index
	for( int m1 = 0 ; m1 < numMoms1 ; m1++)	//time index
		momfile>>mu(m0,m1).real()>>mu(m0,m1).imag();
	momfile.close();

	
	//the momcutoff define the number of moments that are going to be used in the expansion
	//this serve mostly for convergebce test
	int momCutOff = maxNumMom,tmomCutOff = maxNumTMom; //First it is assume that one is going to use all moments
	cfg.lookupValue("MomentsCutOff" ,momCutOff);	//then it will try to find the moments cutoff
	cfg.lookupValue("TMomentCutOff" ,tmomCutOff);	//and the time moments cutoff
	if (momCutOff <= numMoms0 && momCutOff>0 && tmomCutOff <= numMoms1 && tmomCutOff>0 ) //if this are valid number replace the maximum
	{
		maxNumMom= momCutOff; 	maxNumTMom= tmomCutOff;
	}	
	else //If not, trow some warning
		std::cerr<<std::endl<<"Warning."<<std::endl
				 <<"Your MomentsCutOff="<<momCutOff<<std::endl
				 <<"or your TMomentsCutOff="<<tmomCutOff<<" is invalid."<<std::endl
				 <<"Using the maximum values= "<<maxNumMom<<"x"<<maxNumTMom<<" as cutoff"<<std::endl;


	//Read fron config file the energy grid where the conductivity will be computed
	double Emin=-1.0, Emax=1.0, dE=0.1;
	try{
		const libconfig::Setting &enerGrid = cfg.lookup("EnergyGrid");
		std::vector<double> enerGridParam( enerGrid.getLength() );
		for(int i = 0; i < enerGrid.getLength(); ++i)
			enerGridParam[i] = (double)enerGrid [i];
		Emin= enerGridParam[0]; 	
		Emax= enerGridParam[1];
		dE  = enerGridParam[2];
		}
	catch(const libconfig::SettingNotFoundException &nfex)
	{
		std::cerr<<"EnergyGrid field not found, using default values"<<std::endl;
	}
	std::cout	<<"The NonEqRes will be computed in the uniform grid"<<std::endl
				<<"between "<<Emin<<" and "<<Emax<<" with steps "<<dE<<" in units of eV"<<std::endl;

	//Read the temperature, dephasing time, and simulation time 
	const double kb    = 0.086173324; 	// boltlzmann constant meV/K
	const double hbar  = 0.6582119624; 	// reduced  planck constant eV.fs

	//TEMPERATURE
	double temp=300.0; 
	cfg.lookupValue("Temperature",temp); 

	//TIMES ARE HANDLED USING TIME UNIT TO HANDLE INFINITIES
	time_unit tphi = 10.0*temp*kb/hbar ;
	tphi.ReadTime(cfg,"tphi" );
	time_unit tmax; if( !tphi.Infinite() ) tmax = 10.*tphi ; else tmax = 10.0*temp*kb/hbar;
	tmax.ReadTime(cfg,"tmax" );
	if( tphi.Infinite() &&tmax.Infinite() )
	{
		std::cerr<<" tphi and tmax cannot be both infinite simultaneously"<<std::endl;
		return -1; 
	}

	//Scaling  and normalization
	const double tbr   	= temp*kb;		//thermal broadening
	const double ntb   	= tbr/1000./HalfWidth; 	//normalized thermal broadening
	time_unit ntphi 	= tphi*HalfWidth/hbar; //normalized tphi
	time_unit ntmax 	= tmax*HalfWidth/hbar; //normalized tmax
	double eta_tmax = 0.0; if (!tmax.Infinite() ) eta_tmax = 1000*hbar/tmax;//mev
	double eta_phi  = 0.0; if (!tphi.Infinite() ) eta_phi = 1000*hbar/tphi;//mev
	double min_res = eta_phi;
	const double max_inn_sum = (int)( 40.0/(min_res/HalfWidth/1000.) ); 

	std::cout	<<"The calculation will be performed using: "<<std::endl
				<<" temperature = "<<temp<<" K corresponding to "<<hbar/tbr<<" fs or "<<tbr<<" meV"<<std::endl
				<<" tphi = "<<(std::string)tphi<<" fs corresponding to "<<eta_phi<<" meV"<<std::endl
				<<" tmax = "<<(std::string)tmax<<" fs corresponding to "<<eta_tmax<<" meV"<<std::endl
				<<" numberOfMoments = "<<maxNumMom<<" x "<<maxNumTMom<<std::endl
				<<" Minimal energy resolutionis "<<min_res<<" or "<<max_inn_sum<<" chebyshev nodes"<<std::endl<<std::endl; 


	std::string
	outputName  ="NEQ"+S_OPR+LABEL+"KPM_M"+to_string(maxNumMom)+"x"+to_string(maxNumTMom)+"RS"+to_string(numRV)+"tphi"+(std::string)tphi+"tmax"+(std::string)tmax+"c.dat";
	std::cout<<"Saving the data in "<<outputName<<std::endl;

	std::vector<complex> BesselCn_(  maxNumTMom );
	std::ofstream outputfile( outputName.c_str() );


	//One must compute the energy integration weithed by the fermi distribution function
	//for each chemical potential and temperature. In order to explot the 
	//chebyshev polynomial that integration can be done on the chebyshev nodes
	//xn. Therefore, we create an array neq_n which at the end need to be
	// summed. 
	std::vector<double> neq_( maxNumMom,0.0 );

	double lambda_m0 = maxNumMom *eta_phi/HalfWidth/1000.;
	double lambda_mt = maxNumTMom*eta_phi/HalfWidth/1000.;

	//The second step is to goe through the required fermi energies, which in principel should 
	//have not constrains.
	for( double E=Emin; E < Emax; E+=dE)	
	{
		//The fermi energy need to be normalized
		const double nE = (E-BandCenter)/HalfWidth;	
		//then we compute neq_n for all the chebyhsev nodes, this is for energy integration
		for( double n=0.0; n < max_inn_sum; n=n+1.0 )
		{
			//The chebyshev node 
			const double _phin = M_PI*( n + 0.5 )/(double)max_inn_sum;
			const double _xn = cos(_phin); 
			const double _fn = fermi_dirac( ( nE - _xn )/ntb ) ; 
			if( std::fabs(_xn) <= xbound && _fn!= 0.0 )
			{
				#pragma omp parallel 
				{
					const double phin = _phin; 
					const double xn = _xn; 
					const double fn = _fn; 
					//The sum will be performed in parallel. 
					const int 
					tid = omp_get_thread_num(), 
					numThreads = omp_get_num_threads(),
					bsize= (maxNumMom+numThreads-1 ) / numThreads,
					m0o = ( tid   ) * bsize,
					m0f = ( tid +1) * bsize;
					for(int m0=m0o ; m0< m0f; m0++ )
					if( m0 < maxNumMom )
					{ 
						double  my_neq = 0.0;
						for( int mt = 0 ; mt < maxNumTMom ; mt++)
						{
							const complex mom0= ( mu(m0,mt)+std::conj(mu(mt,m0)) )*0.5;
							//kernel
							double g_m0 = 1.0,gn_t=1.0;
							if( true )
							{
								g_m0 = sinh(lambda_m0*( 1.0 - (double)m0/maxNumMom))/sinh(lambda_m0) ;
								gn_t = sinh(lambda_mt*( 1.0 - (double)mt/maxNumTMom))/sinh(lambda_mt) ;
							}
							const double  C_m0 = g_m0*cos( m0*phin )/max_inn_sum;				
							const complex C_nt = gn_t*besselCoeff::dGCoeff( mt, xn ); 
							my_neq += 2.0*fn*C_m0*(C_nt*mom0).imag();
						}
						neq_[m0] += my_neq;
					}
				}
			}
		}
		
		//integrate over normalized energy
		double neq = 0.0 ;
		for( int m0 = 0 ; m0 < maxNumMom ; m0++)
		{
			neq += neq_[m0];
			neq_[m0] = 0.0;
		}
		const double scale_factor = dim/HalfWidth/HalfWidth;
		outputfile<<E<<" "<<neq*scale_factor<<" "<<std::endl;
	}

	return 0;

}
