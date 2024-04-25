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
	const double CUTOFF = 1.00;

class Moments
{
	public:
	typedef std::complex<double>  value_t;
	typedef std::vector< value_t > vector_t;

	//default constructor
	Moments():
	_pNHAM(0),system_label(""),system_size(0),
	band_width(0),band_center(0){};

	//GETTERS
	void getMomentsParams( Moments& mom)
	{
		this->SetHamiltonian( mom.Hamiltonian() ) ; 
		this->SystemLabel( mom.SystemLabel());
		this->BandWidth( mom.BandWidth() );
		this->BandCenter( mom.BandCenter() );
	};	
	
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
	double ScaleFactor() const { return chebyshev::CUTOFF/HalfWidth(); };

	inline
	double ShiftFactor() const { return -BandCenter()/HalfWidth()/chebyshev::CUTOFF; };

	inline 
	vector_t& MomentVector() { return mu ;}

	inline
	value_t& MomentVector(const int i){return  mu[i]; };

	inline
	Moments::vector_t& Chebyshev0(){ return ChebV0; } 

	inline
	Moments::vector_t& Chebyshev1(){ return ChebV1; } 

	inline
	SparseMatrixType& Hamiltonian() const
	{ 
		return *_pNHAM; 
	};


	//SETTERS
	inline
	void SetHamiltonian( SparseMatrixType& NHAM )
	{ 
		if ( this->SystemSize() == 0 ) //Use the rank of the hamiltonian as system size
			this->SystemSize( NHAM.rank() );	
		assert( NHAM.rank() == this->SystemSize()  );
		_pNHAM = &NHAM; 
	};

	//Heavy functions
	int  Rescale2ChebyshevDomain()
	{
		this->Hamiltonian().Rescale(this->ScaleFactor(),this->ShiftFactor());
		return 0;
	};

	inline
	void SetAndRescaleHamiltonian(SparseMatrixType& NHAM)
	{ 
		this->SetHamiltonian(NHAM );
		this->Rescale2ChebyshevDomain();
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

	void SetInitVectors( const vector_t& T0 );

	void SetInitVectors( SparseMatrixType &OP ,const vector_t& T0 );


	inline
	double Rescale2ChebyshevDomain(const double energ)
	{ 
		return (energ - this->BandCenter() )/this->HalfWidth(); 
	};

	int Iterate( );

	//light functions
        int JacksonKernelMomCutOff( const double broad );
	
	//light functions
        double JacksonKernel(const double m,  const double Mom );


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

	Moments1D():_numMoms(0){};

	Moments1D(const size_t m0):_numMoms(m0){ this->MomentVector( Moments::vector_t(_numMoms, 0.0) ); };

	Moments1D( std::string momfilename );

	//GETTERS
	inline
	size_t MomentNumber() const { return _numMoms; };

	inline
	size_t HighestMomentNumber() const { return  _numMoms; };

	//SETTERS

	//OPERATORS
	inline
	Moments::value_t& operator()(const size_t m0) { return this->MomentVector(m0); };

	void MomentNumber(const size_t _numMoms );

	void saveIn(std::string filename);

	//Transformation
	void ApplyJacksonKernel( const double broad );

	// Input/Output
	void Print();

  
	private:
	size_t _numMoms;
};

  
class Moments1D_nonOrth: public Moments1D
{
	public:
  
        Moments1D_nonOrth( const size_t m0 ):Moments1D(m0), S_(NULL){};
  
        Moments1D_nonOrth( SparseMatrixType &S ):Moments1D(), S_(&S){};

        Moments1D_nonOrth(const size_t m0,  SparseMatrixType &S ): Moments1D(m0), S_(&S){};

        Moments1D_nonOrth( std::string& momfilename,  SparseMatrixType &S ): Moments1D(momfilename), S_(&S){};


        void set_S(SparseMatrixType &S) {S_ = &S;};
  
	void SetInitVectors_nonOrthogonal( vector_t& T0 );

	void SetInitVectors_nonOrthogonal( SparseMatrixType &OP ,vector_t& T0 );
  
        int Iterate_nonOrthogonal();

  
	private:

        SparseMatrixType *S_;

};


class Moments2D: public Moments
{
	public:
	Moments2D():numMoms({0,0}){};


        Moments2D(const size_t m0,const size_t m1):numMoms({int(m0),int(m1)}){ this->MomentVector( Moments::vector_t(numMoms[1]*numMoms[0], 0.0) );    };

	Moments2D( std::string momfilename );


	Moments2D( Moments2D mom, const size_t m0,const size_t m1 )
	{
		this->getMomentsParams(mom);
		this->numMoms={int(m0),int(m1)};
		this->MomentVector( Moments::vector_t(numMoms[1]*numMoms[0], 0.0) );    
	};

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

	
	void AddSubMatrix( Moments2D& sub , const int mL, const int mR)
	{
		for(int m0=0; m0<sub.HighestMomentNumber(0); m0++)
		for(int m1=0; m1<sub.HighestMomentNumber(1); m1++)
			this->operator()(mL+m0,mR+m1) += sub(m0,m1);
	} 

	void InsertSubMatrix( Moments2D& sub , const int mL, const int mR)
	{
		for(int m0=0; m0<sub.HighestMomentNumber(0); m0++)
		for(int m1=0; m1<sub.HighestMomentNumber(1); m1++)
			this->operator()(mL+m0,mR+m1) = sub(m0,m1);
	} 
	
	

	void Print();

	private:
	array<int, 2> numMoms;
};

  class MomentsTD : public Moments
  {
	  public:
	  const double HBAR = 0.6582119624 ;//planck constant in eV.fs
	  
	  MomentsTD():
	  _numMoms(1), _maxTimeStep(1),
	  _timeStep(0),
	  _dt(0)
	  {};
	  
	  MomentsTD( const size_t m, const size_t n ): 
	  _numMoms(m), _maxTimeStep(n),
	  _timeStep(0),
	  _dt(0)
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
	  double ChebyshevFreq() const   { return HalfWidth()/chebyshev::CUTOFF/HBAR; };

	  inline
	  double ChebyshevFreq_0() const   { return BandCenter()/HBAR; };
	  
	  //SETTERS
	  
	  void MomentNumber(const size_t mom);
	  
	  void MaxTimeStep(const  size_t maxTimeStep )  {  _maxTimeStep = maxTimeStep; };

	  inline
	  void IncreaseTimeStep(){ _timeStep++; };

	  inline
	  void ResetTime(){ _timeStep=0; };
	  
	  inline
	  void TimeDiff(const double dt ) { _dt = dt; };
	  
	  
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
    double _dt;
  };



  
class Vectors_sliced : public Moments
{
	public: 
	typedef VectorList< Moments::value_t > vectorList_t;


	int NumberOfVectors() const
	{ return numVecs_;}

	int NumberOfSections() const
	{ return num_sections_;}
  
  
        Moments::vector_t& OPV()
        {return OPV_;};
  
  	int SectionSize() const
	{ return section_size_;}

        int LastSectionSize() const
	{ return last_section_size_;}


  	void SetNumberOfVectors( const int x)
	{
	  numVecs_ = x;
	}
  
	int SetNumSections( const int x)
	{

	    num_sections_ = x;
	    section_size_ = SystemSize()/num_sections_; // All of this assuming */* is larger than *%*. Otherwise, i'd change num_sections until it is.
	    last_section_size_ = SystemSize()/num_sections_ + SystemSize() % num_sections_;
	  
	    return 0;
	}
  
        void SetChebmu(const int num_vecs, const size_t last_section_size)
        {
	     Chebmu_(num_vecs, last_section_size);
        }

        vectorList_t& Chebmu(){return Chebmu_;};

	Vectors_sliced()
        {
                SetChebmu(0,0);
                SetNumSections(1);
	};
  

        Vectors_sliced( chebyshev::Vectors_sliced& vec){ *this = vec; };
  
        Vectors_sliced( const int numVecs, const size_t num_sections): numVecs_(numVecs), num_sections_(num_sections), Chebmu_(0,0) {};
  
        Vectors_sliced(Moments& mom, const int numVecs, const size_t num_sections): numVecs_(numVecs), num_sections_(num_sections), Chebmu_(numVecs, mom.SystemSize()/num_sections + mom.SystemSize() % num_sections ) { 

	  section_size_ = mom.SystemSize()/num_sections; // All of this assuming */* is larger than *%*. Otherwise, i'd change num_sections until it is.
	  last_section_size_ = mom.SystemSize()/num_sections + mom.SystemSize() % num_sections;
	  
	  this->getMomentsParams(mom);	  
        };

        Vectors_sliced(Moments2D& mom, const int numVecs, const size_t num_sections): numVecs_(numVecs), num_sections_(num_sections), Chebmu_(numVecs, mom.SystemSize()/num_sections + mom.SystemSize() % num_sections ) { 

	  section_size_ = mom.SystemSize()/num_sections; // All of this assuming */* is larger than *%*. Otherwise, i'd change num_sections until it is.
	  last_section_size_ = mom.SystemSize()/num_sections + mom.SystemSize() % num_sections;
	  
	  this->getMomentsParams(mom);	  
        };

  
  
        int CreateVectorSet()
        {
		try
		{

		        this->SystemSize( this->Hamiltonian().rank() );
			section_size_ = SystemSize()/num_sections_; // All of this assuming */* is larger than *%*. Otherwise, i'd change num_sections until it is.
	                last_section_size_ = SystemSize()/num_sections_ + SystemSize() % num_sections_;

		        const int vec_size  = last_section_size_;
			const int list_size = this->NumberOfVectors();
			Chebmu_.Resize(list_size, vec_size );
		}
		catch (...)
		{ std::cerr<<"Failed to initilize the vector list."<<std::endl;}

		
		return 0;
	}
    
    
	void getMomentsParams( Moments& mom)
	{
		this->SetHamiltonian( mom.Hamiltonian() ) ; 
		this->SystemLabel( mom.SystemLabel());
		this->BandWidth( mom.BandWidth() );
		this->BandCenter( mom.BandCenter() );
		if ( mom.Chebyshev0().size() == this->SystemSize() )
			this->SetInitVectors( mom.Chebyshev0() );
	}

	inline
	size_t Size() const
	{
		return  (long unsigned int)section_size_*
				(long unsigned int)this->NumberOfVectors();
	}

	inline 
	double SizeInGB() const
	{
		const double list_size = this->NumberOfVectors();
		return  sizeof(value_t)*section_size_*list_size*pow(2.0,-30.0);
	}

	inline
	size_t HighestMomentNumber() const { return  this->Chebmu_.ListSize(); };

	inline
	vectorList_t& List() { return this->Chebmu_; };

	inline
	Moments::vector_t& Vector(const size_t m0) { return this->Chebmu_.ListElem(m0); };

	inline
	Moments::value_t& operator()(const size_t m0) { return this->Chebmu_(m0,0); };

        void operator=( Vectors_sliced& vec)
	{
		this->SetHamiltonian( vec.Hamiltonian() ) ; 
		this->SystemLabel( vec.SystemLabel());
		this->BandWidth( vec.BandWidth() );
		this->BandCenter( vec.BandCenter() );
		if ( vec.Chebyshev0().size() == this->SystemSize() )
			this->SetInitVectors( vec.Chebyshev0() );

		num_sections_ = vec.NumberOfSections();
		numVecs_ = vec.NumberOfVectors();
		section_size_ = vec.SectionSize();
	        last_section_size_ = vec.LastSectionSize();

	}

  
	int IterateAllSliced( int );

        int MultiplySliced( SparseMatrixType & , int );

	double MemoryConsumptionInGB();



	private:
	Moments::vector_t OPV_;
        int numVecs_, num_sections_;
        size_t section_size_, last_section_size_;
        vectorList_t Chebmu_;	
  };



class Vectors : public Vectors_sliced
{
	public: 
	typedef VectorList< Moments::value_t > vectorList_t;

  
	
        Vectors(const size_t nMoms,const size_t dim )
        {
	        SetChebmu(nMoms,dim);
	        SetNumSections(1);
	        SetNumberOfVectors(nMoms);
	};

	Vectors( Moments1D& mom )
	{
	        this->getMomentsParams(mom);
		SetNumSections(1);
		SetNumberOfVectors(mom.HighestMomentNumber());
		CreateVectorSet();		
	};
	
        Vectors( Moments2D& mom ):Vectors_sliced(mom, mom.HighestMomentNumber(),  1)
	{
	  	this->getMomentsParams(mom);	        
		SetNumSections(1);
		SetNumberOfVectors(mom.HighestMomentNumber());
		CreateVectorSet();		
	};

	
	Vectors( Moments2D& mom, size_t i  )
	{
		this->getMomentsParams(mom);
		SetNumSections(1);
		SetNumberOfVectors(mom.HighestMomentNumber());
		CreateVectorSet();		
	};    
    
    
        Vectors( MomentsTD& mom )
	{
		this->getMomentsParams(mom);
		SetNumSections(1);
		SetNumberOfVectors(mom.HighestMomentNumber());
		CreateVectorSet();		
	};

  //int IterateAll( );
  //  int Multiply(SparseMatrixType &OP);
  
        int IterateAll( )
        {
	  IterateAllSliced(0);
	  return 0;
	  };
  
	int EvolveAll( const double DeltaT, const double Omega0);
  
        int Multiply( SparseMatrixType &OP )
        {
	  MultiplySliced(OP, 0);
	  return 0;
	  };
  

	double MemoryConsumptionInGB();





};


  
};

#endif 


