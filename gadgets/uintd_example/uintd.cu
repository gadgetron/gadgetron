#include "uintd.h"

template <int LENGTH> __host__ unsigned int& uintd<LENGTH>::operator[] (int i)
{
  return d[i];
}

template struct uintd<1>;
template struct uintd<2>;
template struct uintd<3>;
template struct uintd<4>;
