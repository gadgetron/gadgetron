/*
 * FFT.cpp
 *
 *  Created on: Nov 29, 2011
 *      Author: hansenms
 */

#include "FFT.h"


template <class T> FFT<T>* FFT<T>::instance_ = NULL;

template class FFT<double>;
template class FFT<float>;
