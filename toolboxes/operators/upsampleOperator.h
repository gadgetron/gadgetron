/** \file UpsamplingOperator.h
    \brief Base class for the upsampling operators.

    For instantiation we refer to
    - the class(/file) cuUpsamplingOperator(/.h) for a gpu instantiated operator using the cuNDArray class
*/

#pragma once

#include "linearOperator.h"
#include "vector_td.h"

namespace Gadgetron{
  
  template <class ARRAY_TYPE, unsigned int D> class upsampleOperator
    : public linearOperator<ARRAY_TYPE>
  {
    
  public:
    
    upsampleOperator() : linearOperator<ARRAY_TYPE>() {}
    virtual ~upsampleOperator() {}
    
    typedef typename ARRAY_TYPE::element_type T;

    virtual void mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( accumulate ){
        auto tmp = upsample<T,D>(*in);
        *out += tmp;
      }
      else
        upsample<T,D>(*in,*out);
    }
    
    virtual void mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( accumulate ){
        auto tmp = downsample<T,D>(*in);
        *out += tmp;
      }
      else
        downsample<T,D>(*in,*out);
    }

  };
}
