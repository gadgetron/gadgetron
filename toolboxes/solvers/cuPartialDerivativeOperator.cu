#include "cuPartialDerivativeOperator.h"
#include "vector_td_utilities.h"
#include "check_CUDA.h"
#include "ndarray_vector_td_utilities.h"

template<class T, unsigned int D> __global__ void
first_order_partial_derivative_kernel( typename intd<D>::Type stride, typename intd<D>::Type dims, T *in, T *out )
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
second_order_partial_derivative_kernel( typename intd<D>::Type forwards_stride, typename intd<D>::Type adjoint_stride, typename intd<D>::Type dims, T *in, T *out )
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

template< class REAL, class T, unsigned int D> int 
cuPartialDerivativeOperator<REAL,T,D>::compute_partial_derivative( typename intd<D>::Type stride, cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate )
{
  if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
    std::cerr << std::endl << "partialDerivativeOperator::compute_partial_derivative : array dimensions mismatch." << std::endl;
    return -1;
  }
  if (!accumulate) cuNDA_clear(out);
  
  typename uintd<D>::Type _dims = vector_to_uintd<D>( *(in->get_dimensions().get()) );
  typename intd<D>::Type dims;
  for( unsigned int i=0; i<D; i++ ){
    dims.vec[i] = (int)_dims.vec[i];
  }  
  
  if( D<2 ){
    std::cerr << std::endl << "partialDerivativeOperator::compute_partial_derivative : internal error (only D>1 supported for now)." << std::endl;
    return -1;
  }
  
  _set_device();
  
  dim3 dimBlock( dims.vec[0] );
  dim3 dimGrid( 1, dims.vec[D-1] );
  
  for( unsigned int d=1; d<D-1; d++ )
    dimGrid.x *= dims.vec[d];
  
  // Invoke kernel
  first_order_partial_derivative_kernel<T,D><<< dimGrid, dimBlock >>> ( stride, dims, in->get_data_ptr(), out->get_data_ptr() );
  
  CHECK_FOR_CUDA_ERROR();

  _restore_device();
  
  return 0;
}

template< class REAL, class T, unsigned int D> int 
cuPartialDerivativeOperator<REAL,T,D>::compute_second_order_partial_derivative( typename intd<D>::Type forwards_stride, typename intd<D>::Type adjoint_stride, 
								      cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate )
{
  
  if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
    std::cerr << std::endl << "partialDerivativeOperator::compute_second_order_partial_derivative : array dimensions mismatch." << std::endl;
    return -1;
  }
  if (!accumulate) cuNDA_clear(out);
  typename uintd<D>::Type _dims = vector_to_uintd<D>( *(in->get_dimensions().get()) );
  typename intd<D>::Type dims;
  for( unsigned int i=0; i<D; i++ ){
    dims.vec[i] = (int)_dims.vec[i];
  }  
  
  if( D<2 ){
    std::cerr << std::endl << "partialDerivativeOperator::compute_second_order_partial_derivative : internal error (only D>1 supported for now)." << std::endl;
    return -1;
  }
  
  _set_device();

  dim3 dimBlock( dims.vec[0] );
  dim3 dimGrid( 1, dims.vec[D-1] );
  
  for( unsigned int d=1; d<D-1; d++ )
    dimGrid.x *= dims.vec[d];
  
  // Invoke kernel
  second_order_partial_derivative_kernel<T,D><<< dimGrid, dimBlock >>> ( forwards_stride, adjoint_stride, dims, in->get_data_ptr(), out->get_data_ptr() );
  
  CHECK_FOR_CUDA_ERROR();

  _restore_device();
  
  return 0;
}

// Instantiations

template class EXPORTSOLVERS cuPartialDerivativeOperator<float, float, 2>;
template class EXPORTSOLVERS cuPartialDerivativeOperator<float, float, 3>;
template class EXPORTSOLVERS cuPartialDerivativeOperator<float, float, 4>;

template class EXPORTSOLVERS cuPartialDerivativeOperator<float, float_complext, 2>;
template class EXPORTSOLVERS cuPartialDerivativeOperator<float, float_complext, 3>;
template class EXPORTSOLVERS cuPartialDerivativeOperator<float, float_complext, 4>;

template class EXPORTSOLVERS cuPartialDerivativeOperator<double, double, 2>;
template class EXPORTSOLVERS cuPartialDerivativeOperator<double, double, 3>;
template class EXPORTSOLVERS cuPartialDerivativeOperator<double, double, 4>;

template class EXPORTSOLVERS cuPartialDerivativeOperator<double, double_complext, 2>;
template class EXPORTSOLVERS cuPartialDerivativeOperator<double, double_complext, 3>;
template class EXPORTSOLVERS cuPartialDerivativeOperator<double, double_complext, 4>;
