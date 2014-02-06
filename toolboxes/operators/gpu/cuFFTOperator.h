/** \file cuFFTOperator.h
    \brief Instantiation of the Cartesian FFT operator on the gpu.
    
    The file cuFFTOperator.h is a convienience wrapper for the device independent FFTOperator class.
    The class cuFFTOperator instantiates the FFTOperator for cuNDArray< complext<T> >
    and the header furthermore includes additional neccessary header files.
*/

#pragma once

#include "cuNDArray_math.h"
#include "FFTOperator.h"
#include "cuNDFFT.h"

namespace Gadgetron{
  
  /** \class cuFFTOperator
      \brief Instantiation of the Cartesian FFT operator on the gpu.
      
      The class cuFFTOperator is a convienience wrapper for the device independent FFTOperator.
      It instantiates the FFTOperator for type cuNDArray<T>.
  */
  template <class T> class cuFFTOperator : public FFTOperator< cuNDArray< complext<T> >, cuNDFFT<T> >
  {
  public:    
    cuFFTOperator() : FFTOperator< cuNDArray< complext<T> >, cuNDFFT<T> >() {}
    virtual ~cuFFTOperator() {}
  }; 
}
