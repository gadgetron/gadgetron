/*
 * FFT.cpp
 *
 *  Created on: Nov 29, 2011
 *      Author: hansenms
 */

#include "FFT.h"
namespace Gadgetron{
template <typename T> FFT<T>* FFT<T>::instance()
{
		if (!instance_) instance_ = new FFT<T>();
		return instance_;
}

template <class T> FFT<T>* FFT<T>::instance_ = NULL;

template class FFT<double>;
template class FFT<float>;
}
