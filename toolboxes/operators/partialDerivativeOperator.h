/** \file partialDerivativeOperator.h
    \brief Base class for the partialDerivative operators.

    The file partialDerivativeOperator.h is a device independent partial implementation 
    of a partial derivative operator.
    To simplify the actual instantiation we refer to 
    - the class(/file) hoPartialDerivativeOperator(/.h) for a cpu instantiated operator using the hoNDArray class
    - the class(/file) cuPartialDerivativeOperator(/.h) for a gpu instantiated operator using the cuNDArray class
*/
#pragma once

#include "linearOperator.h"
#include "vector_td.h"

namespace Gadgetron{
  
  template < unsigned int D, class ARRAY_TYPE> class partialDerivativeOperator 
    : public linearOperator<ARRAY_TYPE>
  {
    
  public:
    
    partialDerivativeOperator( size_t dimension ) : 
      linearOperator<ARRAY_TYPE>() { compute_stride(dimension); }
    
    virtual ~partialDerivativeOperator() {}
    
    virtual void mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      compute_partial_derivative( forwards_stride_, in, out, accumulate );
    }
    
    virtual void mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {
      compute_partial_derivative( adjoint_stride_, in, out, accumulate );
    }
    
    virtual void mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
    {    
      compute_second_order_partial_derivative( forwards_stride_, adjoint_stride_, in, out, accumulate );
    }
    
    virtual void compute_partial_derivative
    ( typename int64d<D>::Type stride, 
      ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate ) = 0;  
    
    virtual void compute_second_order_partial_derivative
    ( typename int64d<D>::Type forwards_stride, typename int64d<D>::Type adjoint_stride, 
      ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate ) = 0;    
    
  protected:
    
    virtual void compute_stride( size_t _dimension )
    {
      size_t dim = _dimension;
      
      if( _dimension > D-1 ){
        throw std::runtime_error("Error: partialDerivativeOperator: dimension out of range");
      }
      
      for( unsigned int d=0; d<D; d++ ){
        forwards_stride_.vec[d] = (d==dim) ? 1 : 0;
        adjoint_stride_.vec[d] = (d==dim) ? -1 : 0;
      }    
    }
    
  private:
    typename int64d<D>::Type forwards_stride_;
    typename int64d<D>::Type adjoint_stride_;
  };
}
