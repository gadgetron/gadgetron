/** \file hoPartialDerivativeOperator.h
    \brief Implementation of the partial derivative operator for the cpu.
*/

#pragma once

#include "partialDerivativeOperator.h"
#include "hoNDArray_operators.h"
#include "hoNDArray_elemwise.h"
#include "vector_td_utilities.h"

namespace Gadgetron{

  /** \class hoPartialDerivativeOperator
      \brief CPU implementation of device dependent portions of the partialDerivative operator.
  */
  template <class T, unsigned int D> class hoPartialDerivativeOperator
    : public partialDerivativeOperator<D, hoNDArray<T> >
  {
  public:
    
    hoPartialDerivativeOperator() : 
      partialDerivativeOperator< D, hoNDArray<T> >(0) {}
    
    hoPartialDerivativeOperator( unsigned int dimension ) : 
      partialDerivativeOperator<D, hoNDArray<T> >( dimension ) {}
    
    virtual ~hoPartialDerivativeOperator() {}
    
    virtual void compute_partial_derivative( typename intd<D>::Type stride, hoNDArray<T> *in,
					     hoNDArray<T> *out, bool accumulate )
    {
      if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
	BOOST_THROW_EXCEPTION(runtime_error( "hoPartialDerivativeOperator::compute_partial_derivative : array dimensions mismatch."));	
      }
      
      if( in->get_number_of_dimensions() != D || out->get_number_of_dimensions() != D ){
	BOOST_THROW_EXCEPTION(runtime_error("hoPartialDerivativeOperator::compute_partial_derivative : dimensionality mismatch"));	
      }
      
      typename uintd<D>::Type _dims = vector_to_uintd<D>( *(in->get_dimensions().get()) );
      typename intd<D>::Type dims;
      for( unsigned int i=0; i<D; i++ ){
	dims.vec[i] = (int)_dims.vec[i];
      }  
      
      for( unsigned int idx=0; idx<in->get_number_of_elements(); idx++ ) {
	
	T valN, valC;
	
	typename intd<D>::Type co = idx_to_co<D>(idx, dims);
	typename intd<D>::Type coN = (co+dims+stride)%dims;
	
	valN = in->get_data_ptr()[co_to_idx<D>(coN, dims)];
	valC = in->get_data_ptr()[co_to_idx<D>(co, dims)];
	
	T val = valN-valC;
	
	if( accumulate )
	  out->get_data_ptr()[idx] += val;
	else
	  out->get_data_ptr()[idx] = val;
      }
    }
        
    virtual void compute_second_order_partial_derivative( typename intd<D>::Type forwards_stride,
							  typename intd<D>::Type adjoint_stride, 
							  hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate )
    {
      if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
	BOOST_THROW_EXCEPTION(runtime_error( "hoPartialDerivativeOperator::compute_second_order_partial_derivative : array dimensions mismatch."));	
      }
      
      if( in->get_number_of_dimensions() != D || out->get_number_of_dimensions() != D ){
	BOOST_THROW_EXCEPTION(runtime_error( "hoPartialDerivativeOperator::compute_second_order_partial_derivative : dimensionality mismatch"));	
      }
      
      typename uintd<D>::Type _dims = vector_to_uintd<D>( *(in->get_dimensions().get()) );
      typename intd<D>::Type dims;
      for( unsigned int i=0; i<D; i++ ){
	dims.vec[i] = (int)_dims.vec[i];
      }  
  
      for( unsigned int idx=0; idx<in->get_number_of_elements(); idx++ ) {
	
	T valN1, valN2, valC;
	
	typename intd<D>::Type co = idx_to_co<D>(idx, dims);
	typename intd<D>::Type coN1 = (co+dims+forwards_stride)%dims;
	typename intd<D>::Type coN2 = (co+dims+adjoint_stride)%dims;
	
	valN1 = in->get_data_ptr()[co_to_idx<D>(coN1, dims)];
	valN2 = in->get_data_ptr()[co_to_idx<D>(coN2, dims)];
	valC = in->get_data_ptr()[co_to_idx<D>(co, dims)];
	
	T val = valC+valC-valN1-valN2;
	
	if( accumulate )
	  out->get_data_ptr()[idx] += val;
	else
	  out->get_data_ptr()[idx] = val;
      }
    }

    virtual boost::shared_ptr< linearOperator< hoNDArray<T> > > clone() {
      return linearOperator< hoNDArray<T> >::clone(this);
    }    
  };
}
