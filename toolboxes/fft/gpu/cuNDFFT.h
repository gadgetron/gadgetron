/** \file cuNDFFT.h
    \brief Wrapper of the CUFFT library for ndarrays of type Gadgetron::complext.
*/

#ifndef CUNDFFT_H
#define CUNDFFT_H
#pragma once

#include "cuNDArray.h"

namespace Gadgetron{

  /** \class cuNDFFT
      \brief Wrapper of the CUFFT library for ndarrays of type complext.

      Wrapper of the CUFFT library for ndarrays of type complext<REAL>.
      The class' template type is a REAL, ie. float or double.
      The FFTs are performed in-place.
  */
  template<class T> class cuNDFFT
  {
  public:

    static cuNDFFT<T>* instance();

    void fft ( cuNDArray<complext<T> > *image, std::vector<size_t> *dims_to_transform, bool do_scale = true );
    void ifft( cuNDArray<complext<T> > *image, std::vector<size_t> *dims_to_transform, bool do_scale = true );

    void fft ( cuNDArray<complext<T> > *image, unsigned int dim_to_transform, bool do_scale = true);
    void ifft( cuNDArray<complext<T> > *image, unsigned int dim_to_transform, bool do_scale = true );

    void fft ( cuNDArray<complext<T> > *image, bool do_scale = true );
    void ifft( cuNDArray<complext<T> > *image, bool do_scale = true );

    void fft1(cuNDArray<complext<T> > *image, bool do_scale = true);
    void fft2(cuNDArray<complext<T> > *image, bool do_scale = true);
    void fft3(cuNDArray<complext<T> > *image, bool do_scale = true);

    void ifft1(cuNDArray<complext<T> > *image, bool do_scale = true);
    void ifft2(cuNDArray<complext<T> > *image, bool do_scale = true);
    void ifft3(cuNDArray<complext<T> > *image, bool do_scale = true);


  protected:
    cuNDFFT() {}
    virtual ~cuNDFFT() {}
    void fft_int( cuNDArray<complext<T> > *image, std::vector<size_t> *dims_to_transform, int direction, bool do_scale = true );
    void fft1_int( cuNDArray<complext<T> > *image, int direction, bool do_scale = true );
    void fft2_int( cuNDArray<complext<T> > *image, int direction, bool do_scale = true );
    void fft3_int( cuNDArray<complext<T> > *image, int direction, bool do_scale = true );
    static cuNDFFT<T>* __instance;


  };

  template<class T> void timeswitch(cuNDArray<complext<T> >*, int);
  template<class T> void timeswitch1D(cuNDArray<complext<T> >*);
  template<class T> void timeswitch2D(cuNDArray<complext<T> >*);
  template<class T> void timeswitch3D(cuNDArray<complext<T> >*);
}

#endif
