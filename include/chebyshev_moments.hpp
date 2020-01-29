#ifndef CHEB_MOM
#define CHEB_MOM

#include <array>
#include <vector>
#include <complex>
#include <string>
#include <cassert>
#include <iostream>		/* for std::cout mostly */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <vector>    //for std::vector mostly class
#include <numeric>   //for std::accumulate *
#include <algorithm> //for std::max_elem
#include <complex>   //for std::complex
#include <fstream>   //For ofstream
#include <limits>    //For getting machine precision limit
#include "sparse_matrix.hpp"
#include "linear_algebra.hpp"




using namespace std;

namespace chebyshev 
{


class Moments
{
	public:
	typedef std::complex<double>  value_t;
	typedef std::vector< value_t > vector_t;

	//GETTERS
	inline
	int SystemSize() const { return system_size; };

	inline
	string SystemLabel() const { return system_label; };

	inline
	double BandWidth() const { return band_width; };

	inline
	double HalfWidth() const { return BandWidth()/2.0; };

	inline
	double BandCenter() const { return band_center; };

	inline
	double ScaleFactor() const { return 1/HalfWidth(); };

	inline
	double ShiftFactor() const { return -BandCenter()/HalfWidth(); };

	inline 
	vector_t& MomentVector() { return mu ;}

	inline
	value_t& MomentVector(const int i){return  mu[i]; };


	//SETTERS
	inline
	void SystemSize(const int dim)  { system_size = dim; };

	inline
	void SystemLabel(string label)  { system_label = label; };

	inline
	void BandWidth( const double x)  { band_width = x; };

	inline
	void BandCenter(const double x) { band_center = x; };


	inline 
	void MomentVector(const vector_t _mu ) { mu= _mu;}


	//Heavy functions
	int  Rescale2ChebyshevDomain(SparseMatrixType& H)
	{
		H.Rescale(this->ScaleFactor(),this->ShiftFactor());
		return 0;
	};


	//OPERATORS
	inline
	bool operator == (const Moments& other) const  
	{ 
		return true;
/*			this->system_label	== other.system_label &&
			this->system_size	== other.system_size &&
			this->band_width 	== other.band_width && 
			this->band_center	== other.band_center&&
			this->mu 			== other.mu;*/
	};




	private:
	std::string system_label;
	int system_size;
	double band_width,band_center;
	vector_t mu;	
};


class Moments1D: public Moments
{
	public: 
	Moments1D():numMoms(0){};

	Moments1D(const int m0):numMoms(m0){ assert ( m0>0); this->MomentVector( Moments::vector_t(numMoms, 0.0) ); };


	//GETTERS
	inline
	int MomentNumber() const { return numMoms; };

	inline
	int HighestMomentNumber() const { return  numMoms; };

	//SETTERS



	//OPERATORS
	inline
	Moments::value_t& operator()(const int m0) { return this->MomentVector(m0); };

	inline
	bool operator == (const Moments1D& other) const  
	{ 
		return true;/*
			this->system_label	== other.system_label &&
			this->system_size	== other.system_size &&
			this->numMoms  		== other.numMoms  &&
			this->band_width 	== other.band_width && 
			this->band_center	== other.band_center&&
			this->mu 			== other.mu;*/
	};

	private:
	int numMoms;
};


class Moments2D: public Moments
{
	public: 
	Moments2D():numMoms({0,0}){};
	Moments2D(const int m0,const int m1):numMoms({m0,m1}){ assert ( m0>0 &&m1>0); this->MomentVector( Moments::vector_t(numMoms[1]*numMoms[0], 0.0) );    };


	Moments2D( std::string momfilename );

	//GETTERS

	array<int, 2> MomentNumber() const { return numMoms; };

	int HighestMomentNumber(const int i) const { return numMoms[i]; };

	inline
	int HighestMomentNumber() const { return  (numMoms[1] > numMoms[0]) ? numMoms[1] : numMoms[0]; };

	//SETTERS
	inline
	void MomentNumber(const int mom0, const int mom1 )
	{ 
		assert ( mom0<= numMoms[0] && mom1 <= numMoms[1] );
		Moments2D new_mom( mom0, mom1 );
		for( int m0 = 0 ; m0 < mom0 ; m0++)
		for( int m1 = 0 ; m1 < mom1 ; m1++)
			new_mom(m0,m1) = this->operator()(m0,m1); 
		this->numMoms= new_mom.MomentNumber();
		this->MomentVector( new_mom.MomentVector() );
	};

	//OPERATORS
	inline
	Moments::value_t& operator()(const int m0,const int m1) { return this->MomentVector(m0*numMoms[1] + m1 ); };

	inline
	bool operator == (const Moments2D& other) const  
	{ 
		return true;/*
			this->system_label	== other.system_label &&
			this->system_size	== other.system_size &&
			this->numMoms  		== other.numMoms  &&
			this->band_width 	== other.band_width && 
			this->band_center	== other.band_center&&
			this->mu 			== other.mu;*/
	};


	//Transformation
	void ApplyJacksonKernel()
	{
		const double
		phi_J0 = M_PI/(double)(numMoms[0]+1.0),
		phi_J1 = M_PI/(double)(numMoms[1]+1.0);
		
		double g_D_m0,g_D_m1;
		for( int m0 = 0 ; m0 < numMoms[0] ; m0++)
		{
			g_D_m0=( (numMoms[0]-m0+1)*cos( phi_J0*m0 )+ sin(phi_J0*m0)*cos(phi_J0)/sin(phi_J0) )*phi_J0/M_PI;
			for( int m1 = 0 ; m1 < numMoms[1] ; m1++)
			{ 
				g_D_m1=( (numMoms[1]-m1+1)*cos( phi_J1*m1 )+ sin(phi_J1*m1)*cos(phi_J1)/sin(phi_J1) )*phi_J1/M_PI;
				this->operator()(m0,m1) *= g_D_m0*g_D_m1;
			}
		}
	};
	
	void ApplyJacksonKernel( const double b0, const double b1 )
	{
		assert( b0 >0 && b1>0);
		const double eta0   =  2.0*b0/1000/this->BandWidth();
		const double eta1   =  2.0*b1/1000/this->BandWidth();
		
		int maxMom0=  ceil(M_PI/eta0);
		int maxMom1=  ceil(M_PI/eta1);

		if(  maxMom0 > numMoms[0] ) maxMom0 = numMoms[0];
		if(  maxMom1 > numMoms[1] ) maxMom1 = numMoms[1];
		std::cout<<"Kernel reduced the number of moments to "<<maxMom0<<" "<<maxMom1<<std::endl;
		this->MomentNumber( maxMom0,maxMom1 ) ;



		const double
		phi_J0 = M_PI/(double)(numMoms[0]+1.0),
		phi_J1 = M_PI/(double)(numMoms[1]+1.0);
		
		double g_D_m0,g_D_m1;
		for( int m0 = 0 ; m0 < numMoms[0] ; m0++)
		{
			g_D_m0=( (numMoms[0]-m0+1)*cos( phi_J0*m0 )+ sin(phi_J0*m0)*cos(phi_J0)/sin(phi_J0) )*phi_J0/M_PI;
			for( int m1 = 0 ; m1 < numMoms[1] ; m1++)
			{
				g_D_m1=( (numMoms[1]-m1+1)*cos( phi_J1*m1 )+ sin(phi_J1*m1)*cos(phi_J1)/sin(phi_J1) )*phi_J1/M_PI;
				this->operator()(m0,m1) *= g_D_m0*g_D_m1;
			}
		}
	}
	
	//COSTFUL FUNCTIONS
    void saveIn(std::string filename);

	void Print();
		
	private:
	array<int, 2> numMoms;
};


class Vectors : public Moments
{
	public: 
	Vectors():numMoms(0){};
	Vectors(const int nMoms,const int dim ):numMoms(nMoms){ assert ( nMoms>0); this->SystemSize(dim); this->MomentVector( Moments::vector_t(nMoms*this->SystemSize(), 0.0) );    };
	Vectors( Moments1D mom ):numMoms( mom.HighestMomentNumber() ){ this->SystemSize(mom.SystemSize() ); this->MomentVector( Moments::vector_t(numMoms*this->SystemSize(), 0.0) );    };
	Vectors( Moments2D mom ):numMoms( mom.HighestMomentNumber() ){ this->SystemSize(mom.SystemSize() ); this->MomentVector( Moments::vector_t(numMoms*this->SystemSize(), 0.0) );    };
	Vectors( Moments2D mom, const int i ):numMoms( mom.HighestMomentNumber(i) ){ this->SystemSize(mom.SystemSize() ); this->MomentVector( Moments::vector_t(numMoms*this->SystemSize(), 0.0) );    };


	inline
	Moments::value_t& operator()(const int m0) { return this->MomentVector(m0*this->SystemSize() ); };


	inline
	int HighestMomentNumber() const { return  numMoms; };
	
	inline 
	Moments::vector_t& Chebyshev0(){ return ChebV0; } 

	
	void SetInitVectors( SparseMatrixType &NHAM,const Moments::vector_t& T0 )
	{
		assert( NHAM.rank() == this->SystemSize()&& T0.size() == this->SystemSize() );
		ChebV0 = T0; ChebV1 = T0;
		NHAM.Multiply( ChebV0, ChebV1 );
	};


	void SetInitVectors( SparseMatrixType &NHAM, SparseMatrixType &OP ,const Moments::vector_t& T0 )
	{
		assert( OP.rank() == NHAM.rank() && NHAM.rank() == this->SystemSize()&& T0.size() == this->SystemSize() );

		ChebV0 = T0; ChebV1 = T0;
		  OP.Multiply( ChebV1, ChebV0 );
		NHAM.Multiply( ChebV0, ChebV1 );
	};


	int IterateAll( SparseMatrixType &NHAM )
	{
		const int  dim = NHAM.rank();
		assert( dim == this->SystemSize() );

		linalg::copy(dim, &ChebV0[0],&this->operator()(0) ); 
		for(int m=1; m < this->HighestMomentNumber(); m++ )
		{
			linalg::copy(dim, &ChebV1[0],&this->operator()(m) ); 
			NHAM.Multiply(2.0,ChebV1,-1.0,ChebV0);
			ChebV0.swap(ChebV1);
		}
		
		return 0;
	};



	int Multiply( SparseMatrixType &OP )
	{
		assert( OP.rank() == this->SystemSize() );
		const int  dim = OP.rank();
		
		Moments::vector_t OPV( dim );
		for(int m=0; m < this->HighestMomentNumber(); m++ )
		{
			linalg::copy(dim, &this->operator()(m), &OPV[0] ); 
			OP.Multiply(&OPV[0], &this->operator()(m) );
		}
		
		return 0;
	};


	double MemoryConsumptionInGB()
	{
		return sizeof(value_t)*( this->MomentVector().size()/pow(2.0,30.0)+2.0*this->SystemSize()/pow(2.0,30.0) );
	}
	
	
	private:
	int numMoms;
	Moments::vector_t ChebV0,ChebV1;
};

};

#endif 
