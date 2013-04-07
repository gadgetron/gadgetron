#pragma once
#include "vector_td.h"
#include "cuNDArray.h"
#include "thrust/device_vector.h"
namespace Gadgetron{

template<class T, unsigned int D> void vector_fill(cuNDArray< vector_td<T,D> >* data,  vector_td<T,D> val);
template<class T, unsigned int D> void test_abs(cuNDArray< vector_td<T,D> >* data);
template<class T, unsigned int D> thrust::device_vector<T> test_norm(cuNDArray< vector_td<T,D> >* data);
template<class T, unsigned int D> thrust::device_vector<T> test_min(cuNDArray< vector_td<T,D> >* data);

template<class T, unsigned int D> thrust::device_vector<T> test_max(cuNDArray< vector_td<T,D> >* data);
}
