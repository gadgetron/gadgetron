#include "cuCG.h"


template <class T> cuNDArray<T> cuCG<T>::solve(cuNDArray<T>* rhs)
{
  cuNDArray<T> ret;

  return ret;
}

template class cuCG<float>;
template class cuCG<float2>;
