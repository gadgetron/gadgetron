/** \file laplaceOperator.h
    \brief Base class for the Laplacian operator implementations.
*/

#pragma once

#include "linearOperator.h"

namespace Gadgetron{
  
  template <unsigned int D, class ARRAY_TYPE> class laplaceOperator : public linearOperator<ARRAY_TYPE>
  {    
  public:
    
    laplaceOperator( ) : linearOperator<ARRAY_TYPE>() { }
    virtual ~laplaceOperator() {}
    
    virtual void mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      compute_laplace( in, out, accumulate );
    }
  
    virtual void mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      compute_laplace( in, out, accumulate );
    }
    
  protected:
    virtual void compute_laplace( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate ) = 0;
  };
}
