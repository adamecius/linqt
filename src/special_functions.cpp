#include "special_functions.hpp"

double chebyshev::besselJ(const int n, const double x) {
  
  switch (n) {
  case 0:
    return j0(x);
    break;
  case 1:
    return j1(x);
    break;
  default:
    return jn(n, x);
  }

  return -123456;
}
