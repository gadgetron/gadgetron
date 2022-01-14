#pragma once
#include "cuNDArray.h"
#include "complext.h"
#include "vector_td.h"


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPU__)
#define EXPORTGPUDWT __declspec(dllexport)
#else
#define EXPORTGPUDWT __declspec(dllimport)
#endif
#else
#define EXPORTGPUDWT
#endif

namespace Gadgetron{

//template<class T> void dwt(cuNDArray<T> * in_out);

template<class T, unsigned int D, unsigned int WD> EXPORTGPUDWT void DWT1( cuNDArray<T>* in,
		cuNDArray<T>* out, vector_td<typename realType<T>::Type,WD> wavelet, int dim, int shift = 0);

template<class T, unsigned int D, unsigned int WD> EXPORTGPUDWT void IDWT1( cuNDArray<T>* in,
		cuNDArray<T>* out, vector_td<typename realType<T>::Type,WD> wavelet, int dim, int shift = 0);

}
