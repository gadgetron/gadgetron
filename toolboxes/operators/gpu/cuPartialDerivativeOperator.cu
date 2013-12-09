/** \file cuPartialDerivativeOperator.h
    \brief Implementation of the partial derivative operator for the gpu.
*/

#include "cuPartialDerivativeOperator.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "vector_td_utilities.h"
#include "check_CUDA.h"

namespace Gadgetron{

  template<class T, unsigned long long D> __global__ void
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

  template<class T, unsigned long long D> __global__ void
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

  template< class T, unsigned long long D> void
  cuPartialDerivativeOperator<T,D>::compute_partial_derivative( typename intd<D>::Type stride,
								cuNDArray<T> *in, 
								cuNDArray<T> *out, 
								bool accumulate )
  {
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
      throw std::runtime_error( "partialDerivativeOperator::compute_partial_derivative : array dimensions mismatch.");

    }

    if( in->get_number_of_dimensions() != D || out->get_number_of_dimensions() != D ){
      throw std::runtime_error("partialDerivativeOperator::compute_partial_derivative : dimensionality mismatch");

    }

    if (!accumulate) clear(out);
    
    typename intd<D>::Type dims = to_intd( from_std_vector<unsigned long long,D>( *(in->get_dimensions().get()) ));     
    dim3 dimBlock( dims.vec[0] );
    dim3 dimGrid( 1, dims.vec[D-1] );
  
    for( unsigned int d=1; d<D-1; d++ )
      dimGrid.x *= dims.vec[d];
  
    // Invoke kernel
    first_order_partial_derivative_kernel<T,D><<< dimGrid, dimBlock >>> ( stride, dims, in->get_data_ptr(), out->get_data_ptr() );
  
    CHECK_FOR_CUDA_ERROR();
  }

  template<class T, unsigned long long D> void
  cuPartialDerivativeOperator<T,D>::compute_second_order_partial_derivative( typename intd<D>::Type forwards_stride,
									     typename intd<D>::Type adjoint_stride, 
									     cuNDArray<T> *in, cuNDArray<T> *out, 
									     bool accumulate )
  {  
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
      throw std::runtime_error( "partialDerivativeOperator::compute_second_order_partial_derivative : array dimensions mismatch.");

    }

    if( in->get_number_of_dimensions() != D || out->get_number_of_dimensions() != D ){
      throw std::runtime_error( "partialDerivativeOperator::compute_second_order_partial_derivative : dimensionality mismatch");

    }

    if (!accumulate) clear(out);

    typename intd<D>::Type dims = to_intd( from_std_vector<unsigned long long,D>( *(in->get_dimensions().get()) ));
    dim3 dimBlock( dims.vec[0] );
    dim3 dimGrid( 1, dims.vec[D-1] );
  
    for( unsigned int d=1; d<D-1; d++ )
      dimGrid.x *= dims.vec[d];
  
    // Invoke kernel
    second_order_partial_derivative_kernel<T,D><<< dimGrid, dimBlock >>> ( forwards_stride, adjoint_stride, dims, in->get_data_ptr(), out->get_data_ptr() );
  
    CHECK_FOR_CUDA_ERROR();
  }

  //
  // Instantiations
  //

  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float, 1>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float, 2>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float, 3>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float, 4>;

  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float_complext, 1>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float_complext, 2>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float_complext, 3>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float_complext, 4>;

  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double, 1>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double, 2>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double, 3>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double, 4>;

  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double_complext, 1>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double_complext, 2>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double_complext, 3>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double_complext, 4>;
}
