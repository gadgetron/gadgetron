#pragma once
#include "complext.h"
#include "hoCuNDArray.h"
#include "cuNDArray.h"
#include "gpusolvers_export.h"

namespace Gadgetron{

template<class T> void EXPORTGPUSOLVERS solver_non_negativity_filter(cuNDArray<T>* x , cuNDArray<T>* g);


template<class T> void EXPORTGPUSOLVERS updateF(cuNDArray<T>& data, typename realType<T>::Type alpha ,typename realType<T>::Type sigma);

template<class T> void EXPORTGPUSOLVERS updateFgroup(std::vector<cuNDArray<T> >& datas, typename realType<T>::Type alpha ,typename realType<T>::Type sigma);

}
