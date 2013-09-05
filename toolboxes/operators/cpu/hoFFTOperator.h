/** \file hoFFTOperator.h
    \brief Instantiation of the Cartesian FFT operator on the cpu.
    
    The file hoFFTOperator.h is a convienience wrapper for the device independent FFTOperator class.
    The class hoFFTOperator instantiates the FFTOperator for the hoNDArray< std::complex<T> >
    and the header furthermore includes additional neccessary header files.
*/

#pragma once

#include "hoNDArray_math.h"
#include "FFTOperator.h"
#include "hoFFT.h"

namespace Gadgetron{
  
  /** \class hoFFTOperator
      \brief Instantiation of the Cartesian FFT operator on the cpu.
      
      The class hoFFTOperator is a convienience wrapper for the device independent FFTOperator.
      It instantiates the FFTOperator for type hoNDArray< std::complex<T> >.
  */
  template <class T> class hoFFTOperator : public FFTOperator< hoNDArray< std::complex<T> >, hoFFT<T> >
  {
  public:    
    hoFFTOperator() : FFTOperator< hoNDArray< std::complex<T> >, hoFFT<T> >() {}
    virtual ~hoFFTOperator() {}
  }; 
}
