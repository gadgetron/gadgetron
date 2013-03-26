/*
 * hoFFT.cpp
 *
 *  Created on: Nov 29, 2011
 *      Author: hansenms
 */

#include "hoFFT.h"

namespace Gadgetron{

  template <typename T> hoFFT<T>* hoFFT<T>::instance()
  {
    if (!instance_) instance_ = new hoFFT<T>();
    return instance_;
  }
  
  template <class T> hoFFT<T>* hoFFT<T>::instance_ = NULL;
  
  template class hoFFT<double>;
  template class hoFFT<float>;
}
