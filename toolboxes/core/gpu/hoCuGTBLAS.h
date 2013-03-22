#pragma once
#include "hoCuNDArray.h"
#include "cudaDeviceManager.h"

namespace Gadgetron{
template<class T> T dot(hoCuNDArray<T> *x,hoCuNDArray<T> *y,int device );
template<class T> T dot(hoCuNDArray<T> *x,hoCuNDArray<T> *y){ return dot(x,y,cudaDeviceManager::Instance()->getCurrentDevice()); }

template<class T> typename realType<T>::type nrm2( hoCuNDArray<T>* arr, int device);
template<class T> typename realType<T>::type nrm2( hoCuNDArray<T>* arr){return nrm2(arr,cudaDeviceManager::Instance()->getCurrentDevice());}

template<class T> void axpy(T a, hoCuNDArray<T>* x, hoCuNDArray<T>* y, int device);
template<class T> void axpy(T a, hoCuNDArray<T>* x, hoCuNDArray<T>* y){axpy(a,x,y,cudaDeviceManager::Instance()->getCurrentDevice());}

template<class T> void axpy(T a, hoCuNDArray<complext<T> >* x, hoCuNDArray<complext<T> >* y, int device){
	axpy(complext<T>(a),x,y,device);
}
template<class T> void axpy(T a, hoCuNDArray<complext<T> >* x, hoCuNDArray<complext<T> >* y){axpy(a,x,y,cudaDeviceManager::Instance()->getCurrentDevice());}

/**
 * @brief Gets the index of the index of the element with minimum absolute
 * @param x Input data
 * @param device Device on which to perform computation
 * @return index of absolute minimum values
 */
template<class T> int amin(hoCuNDArray<T>* x,int device);
/**
 * @brief Gets the index of the index of the element with minimum absolute
 * @param x Input data
 * @return index of absolute minimum values
 */
template<class T> int amin(hoCuNDArray<T>* x){return amin(x,cudaDeviceManager::Instance()->getCurrentDevice());}

/**
 * @brief Gets the index of the index of the element with maximum absolute
 * @param x Input data
 * @param device Device on which to perform computation
 * @return index of absolute maximum values
 * @details Note that this returns the C-style index and NOT the Fortran index.
 */
template<class T> int amax(hoCuNDArray<T>* x,int device);
/**
 * @brief Gets the index of the index of the element with maximum absolute
 * @param x Input data
 * @return index of absolute maximum values
 * @details Note that this returns the C-style index and NOT the Fortran index.
 */
template<class T> int amax(hoCuNDArray<T>* x){return amax(x,cudaDeviceManager::Instance()->getCurrentDevice());}

template<class T> typename realType<T>::type asum(hoCuNDArray<T>* x,int device);
template<class T> typename realType<T>::type asum(hoCuNDArray<T>* x){return asum(x,cudaDeviceManager::Instance()->getCurrentDevice());}
}
