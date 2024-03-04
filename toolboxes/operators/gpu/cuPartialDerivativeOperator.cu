/** \file cuPartialDerivativeOperator.h
    \brief Implementation of the partial derivative operator for the gpu.
*/

#include "cuPartialDerivativeOperator.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "vector_td_utilities.h"
#include "check_CUDA.h"

namespace Gadgetron{

  template<class T, unsigned int D> __global__ void
  first_order_partial_derivative_kernel( typename intd<D>::Type stride, 
                                         typename intd<D>::Type dims, 
                                         const T  * __restrict__ in, T * __restrict__ out )
  {
    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
    if( idx < prod(dims) ){

      T valN, valC;

      typename intd<D>::Type co = idx_to_co(idx, dims);
      typename intd<D>::Type coN = (co+dims+stride)%dims;
    
      valN = in[co_to_idx(coN, dims)];
      valC = in[co_to_idx(co, dims)];
    
      T val = valN-valC;
    
      out[idx] += val;
    }
  }

  template<class T, unsigned int D> __global__ void
  second_order_partial_derivative_kernel( typename intd<D>::Type forwards_stride, 
                                          typename intd<D>::Type adjoint_stride, 
                                          typename intd<D>::Type dims, 
                                          const T  * __restrict__ in, T * __restrict__ out )
  {
    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
    if( idx < prod(dims) ){

      T valN1, valN2, valC;

      typename intd<D>::Type co = idx_to_co(idx, dims);
      typename intd<D>::Type coN1 = (co+dims+forwards_stride)%dims;
      typename intd<D>::Type coN2 = (co+dims+adjoint_stride)%dims;
    
      valN1 = in[co_to_idx(coN1, dims)];
      valN2 = in[co_to_idx(coN2, dims)];
      valC = in[co_to_idx(co, dims)];
    
      T val = valC+valC-valN1-valN2;
    
      out[idx] += val;
    }
  }

  template< class T, unsigned int D> void
  cuPartialDerivativeOperator<T,D>::compute_partial_derivative( typename int64d<D>::Type stride,
                                                                cuNDArray<T> *in, 
                                                                cuNDArray<T> *out, 
                                                                bool accumulate )
  {
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
      throw std::runtime_error( "partialDerivativeOperator::compute_partial_derivative : array dimensions mismatch.");

    }


    if (!accumulate) clear(out);
    
    typename int64d<D>::Type dims = vector_td<long long,D>( from_std_vector<size_t,D>( *(in->get_dimensions().get()) ));
    dim3 dimBlock( dims.vec[0] );
    dim3 dimGrid( 1, dims.vec[D-1] );
  
    for(int d=1; d<D-1; d++ )
      dimGrid.x *= dims.vec[d];
  
    size_t elements = in->get_number_of_elements();

    // Invoke kernel
    for (size_t i = 0; i < elements/prod(dims); i++)
    	first_order_partial_derivative_kernel<T,D><<< dimGrid, dimBlock >>> 
        ( vector_td<int,D>(stride), vector_td<int,D>(dims),
          in->get_data_ptr()+i*prod(dims), out->get_data_ptr()+i*prod(dims));
  
    CHECK_FOR_CUDA_ERROR();
  }

  template<class T, unsigned int D> void
  cuPartialDerivativeOperator<T,D>::compute_second_order_partial_derivative( typename int64d<D>::Type forwards_stride,
                                                                             typename int64d<D>::Type adjoint_stride, 
                                                                             cuNDArray<T> *in, cuNDArray<T> *out, 
                                                                             bool accumulate )
  {  
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
      throw std::runtime_error( "partialDerivativeOperator::compute_second_order_partial_derivative : array dimensions mismatch.");
    }
    
    if (!accumulate) clear(out);

    typename int64d<D>::Type dims = vector_td<long long,D>( from_std_vector<size_t,D>( *(in->get_dimensions().get()) ));
    dim3 dimBlock( dims.vec[0] );
    dim3 dimGrid( 1, dims.vec[D-1] );
  
    for(int d=1; d<D-1; d++ )
      dimGrid.x *= dims.vec[d];
  
    size_t elements = in->get_number_of_elements();

    // Invoke kernel
		for (size_t i = 0; i < elements/prod(dims); i++)
			second_order_partial_derivative_kernel<T,D><<< dimGrid, dimBlock >>> 
        ( vector_td<int,D>(forwards_stride), vector_td<int,D>(adjoint_stride), vector_td<int,D>(dims),
          in->get_data_ptr()+i*prod(dims), out->get_data_ptr()+i*prod(dims) );
    
    CHECK_FOR_CUDA_ERROR();
  }

  //
  // Instantiations
  //

  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float, 1>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float, 2>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float, 3>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float, 4>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float, 5>;

  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float_complext, 1>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float_complext, 2>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float_complext, 3>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float_complext, 4>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<float_complext, 5>;

  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double, 1>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double, 2>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double, 3>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double, 4>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double, 5>;

  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double_complext, 1>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double_complext, 2>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double_complext, 3>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double_complext, 4>;
  template class EXPORTGPUOPERATORS cuPartialDerivativeOperator<double_complext, 5>;
}


