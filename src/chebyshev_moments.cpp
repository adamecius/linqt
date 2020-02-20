#include "chebyshev_moments.hpp"

//Heavy functions
int  chebyshev::Moments::Rescale2ChebyshevDomain(SparseMatrixType& H)
{
	H.Rescale(this->ScaleFactor(),this->ShiftFactor());
	return 0;
};
