/** \file UpsamplingOperator.h
    \brief Base class for the upsampling operators.

    For instantiation we refer to
    - the class(/file) cuUpsamplingOperator(/.h) for a gpu instantiated operator using the cuNDArray class
*/

#pragma once

#include "linearOperator.h"
#include "vector_td.h"

namespace Gadgetron{
  
  template < unsigned int D, class ARRAY_TYPE> class upsamplingOperatorOperator 
    : public linearOperator<ARRAY_TYPE>
  {
    
  public:
    
    upsamplingOperatorOperator() : linearOperator<ARRAY_TYPE>() {}
    virtual ~upsamplingOperatorOperator() {}
    
    virtual void mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( accumulate ){
        boost::shared_ptr<ARRAY_TYPE> tmp = upsample(in);
        *out += *tmp;
      }
      else
        upsample(in,out);
    }
    
    virtual void mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( accumulate ){
        boost::shared_ptr<ARRAY_TYPE> tmp = downsample(in);
        *out += *tmp;
      }
      else
        downsample(in,out);
    }
  };
}
