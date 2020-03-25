#ifndef CHEBYSHEV_SOLVER
#define CHEBYSHEV_SOLVER

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
#include "chebyshev_moments.hpp"
#include "linear_algebra.hpp"
#include <omp.h>
#include <chrono>
#include "quantum_states.hpp"
#include "kpm_noneqop.hpp" //Get Batch function
#include "special_functions.hpp"

namespace chebyshev
{
	typedef std::complex<double> value_t;	
	typedef std::vector<value_t> vector_t ;

	namespace sequential
	{

		int DensityExpansionMoments(vector_t& PhiL,vector_t& PhiR,
									 SparseMatrixType &HAM,
									 SparseMatrixType &OP,
									 chebyshev::Moments1D &chebMoms);

		int CorrelationExpansionMoments(const vector_t& PhiL, const vector_t& PhiR,
										SparseMatrixType &HAM,
										SparseMatrixType &OPL,
										SparseMatrixType &OPR,
										chebyshev::Moments2D &chebMoms);

		int ComputeMomTable( chebyshev::Vectors &chebVL,  chebyshev::Vectors& chebVR ,  vector_t& output);

	};


	namespace parallel
	{
		int CorrelationExpansionMoments( 	const int batchSize,
											const vector_t& PhiR, const vector_t& PhiL,
											SparseMatrixType &HAM,
											SparseMatrixType &OPL,
											SparseMatrixType &OPR,  
											chebyshev::Vectors &chevVecL,
											chebyshev::Vectors &chevVecR,
											chebyshev::Moments2D &chebMoms
										);
																					
		int ComputeMomTable( chebyshev::Vectors &chebVL, chebyshev::Vectors & chebVR ,  vector_t& output);

    int TempDensityExpansionMoments(const int batchSize,
				    const double Omega0,
				    const vector_t& Phi,
				    SparseMatrixType &HAM,
				    SparseMatrixType &OP,
				    chebyshev::Vectors &chebVecR, chebyshev::Vectors &chebVecL,
				    chebyshev::MomentsTD &chebMoms);
    
    int ComputeTDMomTable(chebyshev::Vectors &chebVL, chebyshev::Vectors& chebVR, vector_t& output);
	};

	int CorrelationExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OPL, SparseMatrixType &OPR,  chebyshev::Moments2D &chebMoms, StateType type);

	int TempDensityExpansionMoments(int numStates, SparseMatrixType &HAM, SparseMatrixType &OP, SparseMatrixType &PROJ,  chebyshev::MomentsTD &chebMoms, StateType type );

}; // namespace chebyshev

#endif
