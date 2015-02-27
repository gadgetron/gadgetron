/** \file identityOperator.h
    \brief Device independent implementation of the identity operator.

    The file identityOperator.h is a device independent implementation of the identity operator.
    To simplify the actual instantiation we refer to 
    - the class(/file) hoIdentityOperator(/.h) for a cpu instantiated operator using the hoNDArray class
    - the class(/file) cuIdentityOperator(/.h) for a gpu instantiated operator using the cuNDArray class
*/

#pragma once

#include "linearOperator.h"

namespace Gadgetron{

  template <class ARRAY_TYPE> class identityOperator : public linearOperator<ARRAY_TYPE>
  {
  public:
    
    identityOperator() : linearOperator<ARRAY_TYPE>() {}
    virtual ~identityOperator() {}
    
    virtual void mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      if( in == 0x0 || out == 0x0 ){
	throw std::runtime_error("Error: identityOperator::mult_{M,MH,MHM}: illegal array pointer provided");
      }

      // We will do only the most basic dimensionality checking
      if( in->get_number_of_elements() != out->get_number_of_elements() ){
	throw std::runtime_error("Error: identityOperator: in/out dimensions mismatch");
      }
        
      if( accumulate )
    	*out += *in;
      else 
	*out = *in;           
    }
    
    virtual void mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      return mult_M(in, out, accumulate);
    }
    
    virtual void mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      return mult_M(in, out, accumulate);
    }

  };
}
