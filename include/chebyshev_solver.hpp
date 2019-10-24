
// Used for OPENMP functions
#include <omp.h>

// C & C++ libraries
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/
#include <limits>		//For getting machine precision limit

//SPARSE MATRIX
#include "sparse_matrix.hpp"
#include "sp_linalg.hpp"
#include "linalg.hpp"

using namespace std;

namespace chebyshev
{
    struct Configure
    {
        int maxMemory; //bites
        double shift;
        double scaleFactor;
        vector<int> TableSize;
    }

    Class MomTable
    {
        public:
            MomTable(Configure *conf)
            {
                _size = conf->TableSize; 
                int numElems =accumulate_size.begin(),_size.end(), 1,multiplies<int>()); //multiply elements of size
                _data = vector<complex>(numElems); //create a vector to hold these elements
            };

            void SetSystemSize(const int systSize){ systSize=_systSize;}
            complex& operator()(const int m0){ return _data[m0]; }

            vector<int> Size(){return _size: }

            complex& operator()(const int m0,const int m1)
            {   assert( _size.size() == 2 ); 
                return _data[ m1*_size[0] + m0 ];
            }

            vector<int> SizeInDir(const int dir)
            {
                assert( dir < _size.size() );
                return _size[dir]:
            }
            
            int maxSize(){ return  *std::max_element(_size.begin(),_size.end()); }

            void saveIn( std::string outputfile){
                const double bandWidth  = 2*conf->scaleFactor;
                const double bandCenter = -conf->shift*conf->scaleFactor;

                typedef std::numeric_limits< double > dbl;
                outputfile.precision(dbl::digits10);
                outputfile<<_systSize<<" "<<bandWidth<<"  "<<bandCenter<< std::endl;
                //Print the number of moments for all directions in a line
                for(vector<int>::iterator it = _size.begin();
                                          it!= _size.end();  it++ )
                    outputfile<<*it<<" "<< std::endl;
                outputfile<< std::endl;

                for(vector<complex>::iterator   it = _data.begin();
                                                it!= _data.end();  it++ )
                    outputfile<<(*it).real()<<" "<<(*it).imag()<<std::endl;
                outputfile.close();

        private: 
            vector<complex> _data;
            vector<int> _size;
            Configure* _conf;
            int _systSize;
    };

    void CorrelationExpansionMoments(State::Creator& stateCreator, SparseMatrix& HAM,SparseMatrix& OPL,SparseMatrix& OPR, MomTable& cTable )
    {
        //Configure chebyshev
        const int MAX_SIZE = cTable.maxSize();
        cTable.SetSystemSize( HAM.Dim() ); 
        const int 
        DIM = HAM.Dim(),
        MAX_MEM    = cTable.maxMemory,
        VEC_SIZE   = sizeof(complex)*2*DIM; //The two is because one need two set of vectors
        int batchSize  = ( MAX_MEM  + VEC_SIZE -1)/VEC_SIZE ;
        if( batchSize > maxSize )
            batchSize = maxSize;

        const double scalFactor = cTable.scaleFactor;
        const double shift = cTable.shift;
        //Allocate the memory   
        complex* data= new complex[ 2*batchSize*DIM ];
        complex* JL  = &data[0*batchSize*DIM];
        complex* JR  = &data[1*batchSize*DIM];
        
        //INITIALIZE ITERATION
        while ( stateCreator.setNewState)
        {
            Phi = stateCreator.getNewState().begPointer();

            //Start the chebyshev expansion of the correlations
            OPR.Multiply(1.0, Phi  , 0.0 , JR[0] );
            HAM.Multiply(scalFactor, JR[0], 0.0 , JR[1]);
            axpy(DIM,shift,JR[0], JR[1] );
            for( int m0 = 0 ; m0 < cTable.SizeInDir(0); m0+= batchSize )
            {
                for( int mR = m0+2 ; mR < m0+batchSize; mR++ )
                if( mR < cTable.SizeInDir(0) )
                {
                    copy(DIM,JR[mR-2], JR[mR]);
                    HAM.Multiply(2.0*scaleFactor, JR[mR-1], -1.0 , JR[mR] );
                    axpy(DIM,conf.shift,JR[mR-1], JR[mR] );
                }
                
                copy(DIM,Psi, JL[0]);
                HAM.Multiply(scaleFactor, JL[0], 0.0 , JL[1]);
                for( int m1 = 0 ; m1 < cTable.SizeInDir(1); m1+= batchSize)
                {
                    for( int mL = m1+2 ; mL < m1+batchSize; mL++ )
                    if( mL < cTable.SizeInDir(1) )
                    {
                        copy(DIM,JL[mL-2], JL[mL]);
                        HAM.Multiply(2.0*scaleFactor, JL[mL-1], -1.0 , JL[mL] );
                        axpy(DIM,shift,JL[mL-1], JR[mL] );
                    }
                //Compute the moments
                    for( int mR = 0 ; mR <m0+batchSize; mR++)
                    for( int mL = 0 ; mL <m1+batchSize; mL++)
                    if( mR < cTable.SizeInDir(0) && mL < cTable.SizeInDir(1) )
                        cTable(mL,mR) += dot(DIM ,JR[mL], JL[mR] );
                }
            }
        }
        delete [] data;
    }
}