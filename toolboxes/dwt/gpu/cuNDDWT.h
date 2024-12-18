#pragma once
#include "cuNDArray.h"
#include "complext.h"
#include "vector_td.h"


namespace Gadgetron{

//template<class T> void dwt(cuNDArray<T> * in_out);

template<class T, unsigned int D, unsigned int WD> void DWT1( cuNDArray<T>* in,
		cuNDArray<T>* out, vector_td<typename realType<T>::Type,WD> wavelet, int dim, int shift = 0);

template<class T, unsigned int D, unsigned int WD> void IDWT1( cuNDArray<T>* in,
		cuNDArray<T>* out, vector_td<typename realType<T>::Type,WD> wavelet, int dim, int shift = 0);

}
