// C & C++ libraries
#include <cassert>   //for assert
#include <vector>    //for std::vector mostly class
#include <numeric>   //for std::accumulate *
#include <algorithm> //for std::max_elem
#include <complex>   ///for std::complex
#include <fstream> //For ofstream
#include <limits>  //For getting machine precision limit
#include "sparse_matrix.hpp"
#include "linear_algebra.hpp"
#include <cassert>
using namespace std;

namespace chebyshev
{
    struct Configure
    {
	Configure(): maxMemory(0), shift(0), scaleFactor(0){};
        int maxMemory; //bites
        double shift;
        double scaleFactor;
        vector<int> TableSize;
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
        inline complex<double> &operator()(const int m0, const int m1){ assert(size_.size() == 2); return data_[m1 * size_[0] + m0]; };
        inline  int Size_InDir(const int dir) const {assert(dir < size_.size()); return size_[dir]; };
        inline double ScaleFactor() const {  return conf_->scaleFactor; };
        inline double EnergyShift() const { return conf_->shift; };
        chebyshev::Configure *Configure() const { return conf_;};

        //COSTFUL FUNCTIONS
        void saveIn(std::string filename);

    private:
        vector<complex<double> > data_;
        vector<int> size_;
        chebyshev::Configure *conf_;
        int systSize_;
    };

    void CorrelationExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR, MomTable &cTable);

}; // namespace chebyshev
