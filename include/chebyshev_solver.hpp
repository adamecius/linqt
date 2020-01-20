// C & C++ libraries
#include <cassert>   //for assert
#include <array>

#include <vector>    //for std::vector mostly class
#include <numeric>   //for std::accumulate *
#include <algorithm> //for std::max_elem
#include <complex>   ///for std::complex
#include <fstream>   //For ofstream
#include <limits>    //For getting machine precision limit
#include "sparse_matrix.hpp"
#include "linear_algebra.hpp"
#include <cassert>
using namespace std;

namespace chebyshev
{
struct Configure
{
    Configure() : maxMemory(0), shift(0), scaleFactor(0){};
    int maxMemory; //bites
    double shift;
    double scaleFactor;
    vector<int> TableSize;
};


class Moments1D
{
	public: 
	Moments1D():numMoms(0){};

	Moments1D(const int m0):numMoms(m0), mu( vector<complex<double> >(numMoms, 0.0) ){};


	//GETTERS
	inline
	int SystemSize() const { return system_size; };

	inline
	int MomentNumber() const { return numMoms; };

	inline
	int HighestMomentNumber() const { return  numMoms; };

	inline
	string SystemLabel() const { return system_label; };

	inline
	double BandWidth() const { return band_width; };

	inline
	double HalfWidth() const { return BandWidth()/2.0; };

	inline
	double BandCenter() const { return band_center; };


	//SETTERS
	inline
	void SystemSize(const int dim)  { system_size = dim; };

	inline
	void SystemLabel(string label)  { system_label = label; };

	inline
	void BandWidth( const double x)  { band_width = x; };

	inline
	void BandCenter(const double x) { band_center = x; };

	//OPERATORS
	inline
	complex<double>& operator()(const int m0) { return mu[ m0 ]; };

	inline
	bool operator == (const Moments1D& other) const  
	{ 
		return true;
			this->system_label	== other.system_label &&
			this->system_size	== other.system_size &&
			this->numMoms  		== other.numMoms  &&
			this->band_width 	== other.band_width && 
			this->band_center	== other.band_center&&
			this->mu 			== other.mu;
	};

	private:
	std::string system_label;
	int system_size;
	int numMoms;
	double band_width,band_center;
	std::vector< std::complex<double> > mu;

};



class Moments2D
{
	public: 
	Moments2D():numMoms({0,0}){};

	Moments2D(const int m0,const int m1):numMoms({m0,m1}), mu( vector<complex<double> >(numMoms[1]*numMoms[0], 0.0) ){};

	Moments2D( std::string momfilename );

	//GETTERS
	inline
	int SystemSize() const { return system_size; };

	array<int, 2> MomentNumber() const { return numMoms; };

	inline
	void MomentNumber(const int mom0, const int mom1 )
	{ 
		assert ( mom0<= numMoms[0] && mom1 <= numMoms[1] );
		Moments2D new_mom( mom0, mom1 );
		for( int m0 = 0 ; m0 < mom0 ; m0++)
		for( int m1 = 0 ; m1 < mom1 ; m1++)
			new_mom(m0,m1) = this->operator()(m0,m1); 
		this->numMoms= new_mom.MomentNumber();
		this->mu 	 = new_mom.MomentVector();
	};

	int HighestMomentNumber(const int i) const { return numMoms[i]; };

	inline
	int HighestMomentNumber() const { return  (numMoms[1] > numMoms[0]) ? numMoms[1] : numMoms[0]; };

	inline
	string SystemLabel() const { return system_label; };

	inline
	double BandWidth() const { return band_width; };

	inline
	double HalfWidth() const { return BandWidth()/2.0; };

	inline
	double BandCenter() const { return band_center; };


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
	std::vector< std::complex<double> > MomentVector() const { return mu ;}

	//OPERATORS
	inline
	complex<double>& operator()(const int m0,const int m1) { return mu[ m0*numMoms[1] + m1 ]; };

	inline
	bool operator == (const Moments2D& other) const  
	{ 
		return true;
			this->system_label	== other.system_label &&
			this->system_size	== other.system_size &&
			this->numMoms  		== other.numMoms  &&
			this->band_width 	== other.band_width && 
			this->band_center	== other.band_center&&
			this->mu 			== other.mu;
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
		const double eta0   =  2.0*b0/1000/band_width;
		const double eta1   =  2.0*b1/1000/band_width;
		
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
	
	private:
	std::string system_label;
	int system_size;
	array<int, 2> numMoms;
	double band_width,band_center;
	std::vector< std::complex<double> > mu;

};


class MomTable
{
public:
    MomTable(chebyshev::Configure &conf)
    {
        conf_ = &conf;
        size_ = conf_->TableSize;
        int numElems = accumulate(size_.begin(), size_.end(), 1, multiplies<int>()); //multiply elements of size_
        data_ = vector<complex<double> >(numElems);                                   //create a vector to hold these elements
    };

    inline vector<int> Size() const { return size_; };
    inline void SetSystemSize(const int systSize) { systSize_ = systSize; };
    inline complex<double> &operator()(const int m0) { return data_[m0]; };
    inline int maxSize() const { return *std::max_element(size_.begin(), size_.end()); }
    inline complex<double> &operator()(const int m0, const int m1)
    {
        assert(size_.size() == 2);
        return data_[m1 * size_[0] + m0];
    };
    inline int Size_InDir(const int dir) const
    {
        assert(dir < size_.size());
        return size_[dir];
    };
    inline double ScaleFactor() const { return conf_->scaleFactor; };
    inline double EnergyShift() const { return conf_->shift; };
    chebyshev::Configure *Configure() const { return conf_; };

    //COSTFUL FUNCTIONS
    void saveIn(std::string filename);

private:
    vector<complex<double> > data_;
    vector<int> size_;
    chebyshev::Configure *conf_;
    int systSize_;
};



//HELP FUNCTIONS TO BE MOVED
bool GetBatchSize(int& batchSize);

namespace sequential
{
	void CorrelationExpansionMoments(const int numMoms0,const int numMoms1,
									 const std::vector< std::complex<double> >& Phi,
									 SparseMatrixType &HAM, 
									 SparseMatrixType &OPL, 
									 SparseMatrixType &OPR, 
									 chebyshev::MomTable &cTable);
	
};


void CorrelationExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR, MomTable &cTable);

void LocalCorrelationExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR, MomTable &cTable);


void DensityMoments(std::vector< std::complex<double> >& PhiL,
					std::vector< std::complex<double> >& PhiR,
					SparseMatrixType &HAM, 
					SparseMatrixType &OPL, 
					const int numMoms,
					const double scalFactor,
					const double shift,
					std::complex<double>* mu);


void Vectors( 					 
						  std::vector< std::complex<double> >& J0,
						  std::vector< std::complex<double> >& J1,
						 SparseMatrixType &HAM,
						 const int numMoms,
						 const double scalFactor,
						 const double shift,
						 std::vector< std::vector< std::complex<double> > >& JCheb);

}; // namespace chebyshev
