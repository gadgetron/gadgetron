#ifndef CUNDARRAY_NEW_PERMUTE_H
#define CUNDARRAY_NEW_PERMUTE_H

#include "cuNDArray.h"
#include "uintd.h"

template <class T> int cuNDArray_new_permute(cuNDArray<T>* in, cuNDArray<T>* out, std::vector<unsigned int> order);

#endif
