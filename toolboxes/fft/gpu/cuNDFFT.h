/** \file cuNDFFT.h
    \brief Wrapper of the CUFFT library for ndarrays of type Gadgetron::complext.
*/

#ifndef CUNDFFT_H
#define CUNDFFT_H
#pragma once

#include "cuNDArray.h"
#include "gpufft_export.h"

namespace Gadgetron{

  /** \class cuNDFFT
      \brief Wrapper of the CUFFT library for ndarrays of type complext.

      Wrapper of the CUFFT library for ndarrays of type complext<REAL>.
      The class' template type is a REAL, ie. float or double.
      The FFTs are performed in-place.
  */
  template<class T> class EXPORTGPUFFT cuNDFFT
  {
  public:

    static cuNDFFT<T>* instance();

    void fft ( cuNDArray<complext<T> > *image, std::vector<size_t> *dims_to_transform );
    void ifft( cuNDArray<complext<T> > *image, std::vector<size_t> *dims_to_transform, bool do_scale = true );

    void fft ( cuNDArray<complext<T> > *image, unsigned int dim_to_transform);
    void ifft( cuNDArray<complext<T> > *image, unsigned int dim_to_transform, bool do_scale = true );

    void fft ( cuNDArray<complext<T> > *image );
    void ifft( cuNDArray<complext<T> > *image, bool do_scale = true );

  protected:   
    cuNDFFT() {}
    virtual ~cuNDFFT() {}
    void fft_int( cuNDArray<complext<T> > *image, std::vector<size_t> *dims_to_transform, int direction, bool do_scale = true );
    static cuNDFFT<T>* __instance;
  };
}

#endif
