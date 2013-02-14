#pragma once

#include "linearOperator.h"
#include "vector_td.h"

namespace Gadgetron{
template < unsigned int D, class ARRAY_TYPE> class partialDerivativeOperator
	: public linearOperator<ARRAY_TYPE>
{
  
public:
  
  partialDerivativeOperator( unsigned int dimension ) : 
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
  ( typename intd<D>::Type stride, ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate ) = 0;  

  virtual void compute_second_order_partial_derivative
  ( typename intd<D>::Type forwards_stride, typename intd<D>::Type adjoint_stride, 
    ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate ) = 0;    

protected:
  
  virtual void compute_stride( unsigned int _dimension )
  {
    unsigned int dim = _dimension;

    if( _dimension > D-1 ){
      std::cerr << std::endl << "Warning: partialDerivativeOperator::compute_stride : dimension out of range, clamping." << std::endl;
      dim = D-1;
    }
    
    for( unsigned int d=0; d<D; d++ ){
      forwards_stride_.vec[d] = (d==dim) ? 1 : 0;
      adjoint_stride_.vec[d] = (d==dim) ? -1 : 0;
    }
    
  }
  
private:
  typename intd<D>::Type forwards_stride_;
  typename intd<D>::Type adjoint_stride_;
};
}
