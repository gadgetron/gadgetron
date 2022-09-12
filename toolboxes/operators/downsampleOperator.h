/** \file DownsamplingOperator.h
    \brief Base class for the downsampling operators.

    For instantiation we refer to
    - the class(/file) cuDownsamplingOperator(/.h) for a gpu instantiated operator using the cuNDArray class
*/

#pragma once

#include "linearOperator.h"
#include "vector_td.h"

namespace Gadgetron{
  
  template <class ARRAY_TYPE, unsigned int D> class downsampleOperator
    : public linearOperator<ARRAY_TYPE>
  {
    
  public:

    typedef typename ARRAY_TYPE::element_type T;

    downsampleOperator() : linearOperator<ARRAY_TYPE>() {}
    virtual ~downsampleOperator() {}
    
    virtual void mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( accumulate ){
        auto tmp = downsample<T,D>(*in);
        *out += tmp;
      }
      else
        downsample<T,D>(*in,*out);
    }
    
    virtual void mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( accumulate ){
        auto tmp = upsample<T,D>(*in);
        *out += tmp;
      }
      else
        upsample<T,D>(*in,*out);
    }

  };
}
