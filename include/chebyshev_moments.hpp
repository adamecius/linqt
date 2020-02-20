#ifndef CHEBYSHEV_MOMENTS
#define CHEBYSHEV_MOMENTS


#include <complex>
#include <vector>
#include <string>
#include <array>

#include "sparse_matrix.hpp" //contain SparseMatrixType
#include <cassert>			 //needed for assert
#include <fstream>   		 //For ifstream and ofstream
#include <limits>    		 //Needed for dbl::digits10
#include "linear_algebra.hpp"

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
	double ScaleFactor() const { return 1.0/HalfWidth(); };

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
	int  Rescale2ChebyshevDomain(SparseMatrixType& H);


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
	void MomentNumber(const int mom0, const int mom1 );

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
	void ApplyJacksonKernel( const double b0, const double b1 );
	
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

	
	void SetInitVectors( SparseMatrixType &NHAM,const Moments::vector_t& T0 );


	void SetInitVectors( SparseMatrixType &NHAM, SparseMatrixType &OP ,const Moments::vector_t& T0 );


	int IterateAll( SparseMatrixType &NHAM );



	int Multiply( SparseMatrixType &OP );


	double MemoryConsumptionInGB();
	
	
	private:
	int numMoms;
	Moments::vector_t ChebV0,ChebV1;
};

};

#endif 



//#include <array>
//#include <cassert>
//#include <iostream>		/* for std::cout mostly */
//#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
//#include <vector>    //for std::vector mostly class
//#include <numeric>   //for std::accumulate *
//#include <algorithm> //for std::max_elem
//#include <complex>   //for std::complex
//#include "mkl.h"
