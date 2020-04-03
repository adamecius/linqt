// With contributions made by Angel D. Prieto S.
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
#include "vector_list.hpp"
#include "special_functions.hpp"
#include "chebyshev_coefficients.hpp"

namespace chebyshev 
{


class Moments
{
	public:
	typedef std::complex<double>  value_t;
	typedef std::vector< value_t > vector_t;

	//GETTERS
	inline
	size_t SystemSize() const { return system_size; };

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

	inline
	Moments::vector_t& Chebyshev0(){ return ChebV0; } 

	inline
	Moments::vector_t& Chebyshev1(){ return ChebV1; } 

	inline
	SparseMatrixType& Hamiltonian()
	{ 
		return *_pNHAM; 
	};


	//SETTERS
	inline
	void SetHamiltonian(SparseMatrixType& NHAM)
	{ 
		this->Rescale2ChebyshevDomain(NHAM);
		_pNHAM = &NHAM; 
	};

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

//	void SetInitVectors( SparseMatrixType &NHAM,const Moments::vector_t& T0 );//to be deprecated 

	void SetInitVectors( SparseMatrixType &NHAM, SparseMatrixType &OP ,const vector_t& T0 );//to be deprecated 

	void SetInitVectors( const vector_t& T0 );

	void SetInitVectors( SparseMatrixType &OP ,const vector_t& T0 );


	inline
	double Rescale2ChebyshevDomain(const double energ)
	{ 
		return (energ - this->BandCenter() )/this->HalfWidth(); 
	};

	int Iterate( SparseMatrixType &NHAM );

	int Iterate( );

	//light functions
    int JacksonKernelMomCutOff( const double broad );
	
	//light functions
    double JacksonKernel(const double m,  const double Mom );

	//Heavy functions
	int  Rescale2ChebyshevDomain(SparseMatrixType& H);

	private:
	SparseMatrixType* _pNHAM;

	Moments::vector_t ChebV0,ChebV1,OPV;
	std::string system_label;
	size_t system_size;
	double band_width,band_center;
	vector_t mu;	
};


class Moments1D: public Moments
{
	public: 

	Moments1D():numMoms(0){};

	Moments1D(const size_t m0):numMoms(m0){ this->MomentVector( Moments::vector_t(numMoms, 0.0) ); };

	//GETTERS
	inline
	size_t MomentNumber() const { return numMoms; };

	inline
	size_t HighestMomentNumber() const { return  numMoms; };

	//SETTERS

	//OPERATORS
	inline
	Moments::value_t& operator()(const size_t m0) { return this->MomentVector(m0); };

	//Transformation
//	void ApplyJacksonKernel( const double b0 );

	// Input/Output
	void Print();

	private:
	size_t numMoms;
};


class Moments2D: public Moments
{
	public:
	Moments2D():numMoms({0,0}){};

	Moments2D(const size_t m0,const size_t m1):numMoms({m0,m1}){ this->MomentVector( Moments::vector_t(numMoms[1]*numMoms[0], 0.0) );    };

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

	//Transformation
	void ApplyJacksonKernel( const double b0, const double b1 );

	// Input/Output 
	//COSTFUL FUNCTIONS
	void saveIn(std::string filename);

	void Print();

	private:
	array<int, 2> numMoms;
};

  class MomentsTD : public Moments
  {
	  public:
	  
	  MomentsTD():
	  _numMoms(1), _maxTimeStep(1),
	  _timeStep(0),
	  _dt(0), _omega(0)
	  {};
	  
	  MomentsTD( const size_t m, const size_t n ): 
	  _numMoms(m), _maxTimeStep(n),
	  _timeStep(0),
	  _dt(0), _omega(0)
	  { this->MomentVector( Moments::vector_t(m*n, 0.0) );    };
	  
	  MomentsTD( std::string momfilename );
	  
	  //GETTERS
	  inline
	  size_t MomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t HighestMomentNumber() const { return _numMoms;};
	  
	  inline
	  size_t CurrentTimeStep() const { return _timeStep;};

	  inline 
	  size_t MaxTimeStep() const { return _maxTimeStep; };
	  
	  inline
	  double TimeDiff() const   { return _dt; };
	  	  
	  inline
	  double ChebyshevFreq() const   { return _omega; };
	  
	  //SETTERS
	  void MomentNumber(const size_t mom);
	  
	  void MaxTimeStep(const  size_t maxTimeStep )  {  _maxTimeStep = maxTimeStep; };

	  inline
	  void IncreaseTimeStep(){ _timeStep++; };
	  
	  inline
	  void TimeDiff(const double dt ) { _dt = dt; };
	  
	  inline
	  void ChebyshevFreq(const double omega) { _omega = omega; };
	  
	  int Evolve(  vector_t& Phi);
	  
	  //OPERATORS
	  inline
	  Moments::value_t& operator()(const size_t m, const size_t n)
	  {
		  return this->MomentVector( m*MaxTimeStep() + n );
	  };
	  
	  //Transformation
	  void ApplyJacksonKernel( const double broad );
	  
	  //COSTFUL FUNCTIONS
	  void saveIn(std::string filename);
	  
	  void Print();

  private:
    size_t _numMoms, _maxTimeStep, _timeStep;
    double _dt, _omega;
  };

class Vectors : public Moments
{
	public: 
	typedef VectorList< Moments::value_t > vectorList_t;

	Vectors():Chebmu(0,0){};
	Vectors(const size_t nMoms,const size_t dim ):Chebmu(nMoms,dim) {  };
	Vectors( Moments1D mom ): Chebmu(mom.HighestMomentNumber(),mom.SystemSize() ) {  };
	Vectors( Moments2D mom ): Chebmu(mom.HighestMomentNumber(),mom.SystemSize() ) {  };
	Vectors( Moments2D mom, const size_t i ): Chebmu(mom.HighestMomentNumber(i), mom.SystemSize() ) {  };
    Vectors( MomentsTD mom ): Chebmu(mom.HighestMomentNumber(),mom.SystemSize() ) {  };

	size_t Size() const
	{
		return  (long unsigned int)this->Chebmu.VectorSize()*
				(long unsigned int)this->Chebmu.ListSize();
	}

	size_t SystemSize() const { return this->Chebmu.VectorSize(); }

	inline
	size_t HighestMomentNumber() const { return  this->Chebmu.ListSize(); };


	inline
	vectorList_t& List() { return this->Chebmu; };

	inline
	Moments::vector_t& Vector(const size_t m0) { return this->Chebmu.ListElem(m0); };

	inline
	Moments::value_t& operator()(const size_t m0) { return this->Chebmu(m0,0); };

	int IterateAll( SparseMatrixType &NHAM );

	int EvolveAll( SparseMatrixType &NHAM, const double DeltaT, const double Omega0);

	int Multiply( SparseMatrixType &OP );


	double MemoryConsumptionInGB();


	private:
	Moments::vector_t ChebV0,ChebV1,OPV;
	vectorList_t Chebmu;	

};

};

#endif 


