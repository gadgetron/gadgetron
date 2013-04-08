#pragma once

#include "cuNDArray_operators.h"
#include "cuNDArray_utils.h"
#include "NFFT_utils.h"
#include "cuNDFFT.h"
#include "vector_td_utilities.h"
#include "convolutionOperator.h"
#include "gpuoperators_export.h"

namespace Gadgetron{

  template <class REAL, unsigned int D> class EXPORTGPUOPERATORS cuConvolutionOperator 
    : public convolutionOperator<cuNDArray<complext<REAL> >, D >
  {
    
  public:
  
    cuConvolutionOperator() : convolutionOperator<cuNDArray<complext<REAL> >, D>() {  }
    virtual ~cuConvolutionOperator() {}
        
    virtual void operator_fft( bool forwards_transform, cuNDArray< complext<REAL> > *image );
    virtual void origin_mirror( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out );
    
    virtual boost::shared_ptr< linearOperator<cuNDArray< complext<REAL> > > > clone()
    {
      return linearOperator< cuNDArray< complext<REAL> > >::clone(this);
    }
  };
}
