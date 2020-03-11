#include "special_functions.hpp"

double chebyshev::besselJ(const int n, const double x) {
  double y;

  switch (n) {
  case 0:
    y = j0(x);
    break;
  case 1:
    y = j1(x);
    break;
  default:
    y = jn(n, x);
  }

  return y;
}
