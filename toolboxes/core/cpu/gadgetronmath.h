#pragma once

/** \file gadgetronmath.h
\brief Math utility functionx

*/

#define _USE_MATH_DEFINES
#include <math.h>

namespace Gadgetron {

template <typename T> T sinc(T x) {
  
  T val;
  if (std::abs(x)<.01) {
    // to 6th order
    val = 1.0 - 1/6.*std::pow(M_PI*x,2) + 1/120.*std::pow(M_PI*x,4) - 1/5040.*std::pow(M_PI*x,6);
  } else {
    val = std::sin(M_PI*x) / (M_PI*x);
  }

  return val;
}

}
