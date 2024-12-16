/** \file cuConvolutionOperator.h
    \brief Convolution operator, GPU based.
*/

#pragma once

#include "cuNDArray_math.h"
#include "cuNDFFT.h"
#include "vector_td_utilities.h"
#include "convolutionOperator.h"

namespace Gadgetron{

  template <class REAL, unsigned int D> class cuConvolutionOperator
    : public convolutionOperator<cuNDArray<complext<REAL> >, D >
  {

  public:

    cuConvolutionOperator() : convolutionOperator<cuNDArray<complext<REAL> >, D>() {  }
    virtual ~cuConvolutionOperator() {}

    virtual void operator_fft( bool forwards_transform, cuNDArray< complext<REAL> > *image );
    virtual void origin_mirror( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out );

  };
}
