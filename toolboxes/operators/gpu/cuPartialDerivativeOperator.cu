#include "cuPartialDerivativeOperator.h"
#include "vector_td_utilities.h"
#include "check_CUDA.h"
#include "ndarray_vector_td_utilities.h"

namespace Gadgetron{
  template<class T, unsigned int D> __global__ void
  first_order_partial_derivative_kernel( typename intd<D>::Type stride, 
					 typename intd<D>::Type dims, 
					 T *in, T *out )
  {
    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
    if( idx < prod(dims) ){

      T valN, valC;

      typename intd<D>::Type co = idx_to_co<D>(idx, dims);
      typename intd<D>::Type coN = (co+dims+stride)%dims;
    
      valN = in[co_to_idx<D>(coN, dims)];
      valC = in[co_to_idx<D>(co, dims)];
    
      T val = valN-valC;
    
      out[idx] += val;
    }
  }

  template<class T, unsigned int D> __global__ void
  second_order_partial_derivative_kernel( typename intd<D>::Type forwards_stride, 
					  typename intd<D>::Type adjoint_stride, 
					  typename intd<D>::Type dims, 
					  T *in, T *out )
  {
    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
    if( idx < prod(dims) ){

      T valN1, valN2, valC;

      typename intd<D>::Type co = idx_to_co<D>(idx, dims);
      typename intd<D>::Type coN1 = (co+dims+forwards_stride)%dims;
      typename intd<D>::Type coN2 = (co+dims+adjoint_stride)%dims;
    
      valN1 = in[co_to_idx<D>(coN1, dims)];
      valN2 = in[co_to_idx<D>(coN2, dims)];
      valC = in[co_to_idx<D>(co, dims)];
    
      T val = valC+valC-valN1-valN2;
    
      out[idx] += val;
    }
  }

  template< class T, unsigned int D> void
  cuPartialDerivativeOperator<T,D>::compute_partial_derivative( typename intd<D>::Type stride,
								cuNDArray<T> *in, 
								cuNDArray<T> *out, 
								bool accumulate )
  {
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
      BOOST_THROW_EXCEPTION(runtime_error( "partialDerivativeOperator::compute_partial_derivative : array dimensions mismatch."));

    }

    if( in->get_number_of_dimensions() != D || out->get_number_of_dimensions() != D ){
      BOOST_THROW_EXCEPTION(runtime_error("partialDerivativeOperator::compute_partial_derivative : dimensionality mismatch"));

    }

    typename uintd<D>::Type _dims = vector_to_uintd<D>( *(in->get_dimensions().get()) );
    typename intd<D>::Type dims;
    for( unsigned int i=0; i<D; i++ ){
      dims.vec[i] = (int)_dims.vec[i];
    }  
    
    _set_device();
  
    if (!accumulate) out->clear();
  
    dim3 dimBlock( dims.vec[0] );
    dim3 dimGrid( 1, dims.vec[D-1] );
  
    for( unsigned int d=1; d<D-1; d++ )
      dimGrid.x *= dims.vec[d];
  
    // Invoke kernel
    first_order_partial_derivative_kernel<T,D><<< dimGrid, dimBlock >>> ( stride, dims, in->get_data_ptr(), out->get_data_ptr() );
  
    CHECK_FOR_CUDA_ERROR();

    _restore_device();
  }

  template<class T, unsigned int D> void
  cuPartialDerivativeOperator<T,D>::compute_second_order_partial_derivative( typename intd<D>::Type forwards_stride,
									     typename intd<D>::Type adjoint_stride, 
									     cuNDArray<T> *in, cuNDArray<T> *out, 
									     bool accumulate )
  {  
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
      BOOST_THROW_EXCEPTION(runtime_error( "partialDerivativeOperator::compute_second_order_partial_derivative : array dimensions mismatch."));

    }

    if( in->get_number_of_dimensions() != D || out->get_number_of_dimensions() != D ){
      BOOST_THROW_EXCEPTION(runtime_error( "partialDerivativeOperator::compute_second_order_partial_derivative : dimensionality mismatch"));

    }

    typename uintd<D>::Type _dims = vector_to_uintd<D>( *(in->get_dimensions().get()) );
    typename intd<D>::Type dims;
    for( unsigned int i=0; i<D; i++ ){
      dims.vec[i] = (int)_dims.vec[i];
    }  
  
    _set_device();

    if (!accumulate) out->clear();

    dim3 dimBlock( dims.vec[0] );
    dim3 dimGrid( 1, dims.vec[D-1] );
  
    for( unsigned int d=1; d<D-1; d++ )
      dimGrid.x *= dims.vec[d];
  
    // Invoke kernel
    second_order_partial_derivative_kernel<T,D><<< dimGrid, dimBlock >>> ( forwards_stride, adjoint_stride, dims, in->get_data_ptr(), out->get_data_ptr() );
  
    CHECK_FOR_CUDA_ERROR();

    _restore_device();
  }

  //
  // Instantiations
  //

  template class EXPORTSOLVERS cuPartialDerivativeOperator<float, 1>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<float, 2>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<float, 3>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<float, 4>;

  template class EXPORTSOLVERS cuPartialDerivativeOperator<float_complext, 1>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<float_complext, 2>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<float_complext, 3>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<float_complext, 4>;

  template class EXPORTSOLVERS cuPartialDerivativeOperator<double, 1>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<double, 2>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<double, 3>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<double, 4>;

  template class EXPORTSOLVERS cuPartialDerivativeOperator<double_complext, 1>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<double_complext, 2>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<double_complext, 3>;
  template class EXPORTSOLVERS cuPartialDerivativeOperator<double_complext, 4>;
}
