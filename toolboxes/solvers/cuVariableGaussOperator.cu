#include "cuVariableGaussOperator.h"
#include "math_constants.h"
#include "check_CUDA.h"
#include "vector_td_utilities.h"

#define BLOCK_SIZE 512
template<class REAL, class T, unsigned int D> __global__ void
mult_M_kernel( typename intd<D>::Type dims, T *in, T *out,REAL *sigma, REAL *norm )
{  
  const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < prod(dims) ){
    
    __shared__ REAL shared[BLOCK_SIZE];
    __shared__ REAL sSigma[BLOCK_SIZE];
    __shared__ REAL sNorm[BLOCK_SIZE];

    T val = T(0);
    //REAL s = 1.0/(2.0*5.0*5.0);
    REAL s;

    typename intd<D>::Type co2;
    typename intd<D>::Type co = idx_to_co<D>(idx, dims);
    
    for (int k = 0; k < gridDim.x; k++){
      shared[threadIdx.x] = in[k*blockDim.x + threadIdx.x];
      sSigma[threadIdx.x] = sigma[k*blockDim.x + threadIdx.x];
      sNorm[threadIdx.x] = norm[k*blockDim.x + threadIdx.x];
      __syncthreads();

      for (int i = 0; i < blockDim.x; i++){
	s = sSigma[i];

	co2 = idx_to_co<D>(k*blockDim.x+i, dims)-co;
	val += shared[i]*sNorm[i]*exp(- ((REAL)dot<int,D>(co2,co2))*s*s*0.5);		  
      }
    }
    out[idx] = val;    
  }   
}

template<class REAL, class T, unsigned int D> __global__ void
norm_kernel( typename intd<D>::Type dims, REAL *sigma, REAL *out )
{
  const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < prod(dims) ){
    
    __shared__ REAL shared[BLOCK_SIZE];

    T val = T(0);
    REAL s = sigma[idx];

    typename intd<D>::Type co2;
    typename intd<D>::Type co = idx_to_co<D>(idx, dims);
    
    for (int k = 0; k < gridDim.x; k++){
      for (int i = 0; i < blockDim.x; i++){
	co2 = idx_to_co<D>(k*blockDim.x+i, dims)-co;
	val += exp(- ((REAL)dot<int,D>(co2,co2))*s*s*0.5);		  
      }
    }
    out[idx] = 1.0/val;    
  }    
}

template<class REAL, class T, unsigned int D> __global__ void
mult_MH_kernel( typename intd<D>::Type dims, T *in, T *out,REAL *sigma,REAL *norm )
{
  const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < prod(dims) ){
    
    __shared__ REAL shared[BLOCK_SIZE];

    T val = T(0);
    REAL s = sigma[idx];

    typename intd<D>::Type co2;
    typename intd<D>::Type co = idx_to_co<D>(idx, dims);
    
    for (int k = 0; k < gridDim.x; k++){
      shared[threadIdx.x] = in[k*blockDim.x + threadIdx.x];
      __syncthreads();

      for (int i = 0; i < blockDim.x; i++){
	co2 = idx_to_co<D>(k*blockDim.x+i, dims)-co;
	val += shared[i]*exp(- ((REAL)dot<int,D>(co2,co2))*s*s*0.5);		  
      }
    }
    out[idx] = val*norm[idx];    
  }    
}

template< class REAL, class T, unsigned int D> void
cuVariableGaussOperator<REAL,T,D>::set_sigma( cuNDArray<REAL> *sigma )
{
  _sigma = sigma;
  
  typename uintd<D>::Type _dims = vector_to_uintd<D>( *(_sigma->get_dimensions().get()) );
  typename intd<D>::Type dims;
  for( unsigned int i=0; i<D; i++ ){
    dims.vec[i] = (int)_dims.vec[i];
  }  

  dim3 dimBlock( BLOCK_SIZE );
  dim3 dimGrid( prod(dims)/BLOCK_SIZE );
  
  _norm = boost::shared_ptr<cuNDArray<REAL> >(new cuNDArray<REAL>);
  _norm->create(_sigma->get_dimensions().get());
  
  // Invoke kernel
  norm_kernel<REAL,T,D><<< dimGrid, dimBlock >>> (dims, sigma->get_data_ptr(), _norm->get_data_ptr() );
  
  CHECK_FOR_CUDA_ERROR();
}

template< class REAL, class T, unsigned int D> int
cuVariableGaussOperator<REAL,T,D>::mult_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate )
{
  if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
    std::cerr << std::endl << "laplaceOperator::compute_laplace : array dimensions mismatch." << std::endl;
    return -1;
  }
  
  typename uintd<D>::Type _dims = vector_to_uintd<D>( *(in->get_dimensions().get()) );
  typename intd<D>::Type dims;
  for( unsigned int i=0; i<D; i++ ){
    dims.vec[i] = (int)_dims.vec[i];
  }  
  
  if( D>3 ){
    std::cerr << std::endl << "partialDerivativeOperator::compute_laplace : internal error (only D<4 supported for now)." << std::endl;
    return -1;
  }

  _set_device();
  
  dim3 dimBlock( BLOCK_SIZE );
  dim3 dimGrid( prod(dims)/BLOCK_SIZE );
    
  // Invoke kernel
  mult_M_kernel<REAL,T,D><<< dimGrid, dimBlock >>> (dims, in->get_data_ptr(), out->get_data_ptr(), _sigma->get_data_ptr() ,_norm->get_data_ptr() );
  
  CHECK_FOR_CUDA_ERROR();

  _restore_device();

  return 0;
}

template< class REAL, class T, unsigned int D> int
cuVariableGaussOperator<REAL,T,D>::mult_MH( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate)
{
  if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
    std::cerr << std::endl << "laplaceOperator::compute_laplace : array dimensions mismatch." << std::endl;
    return -1;
  }
  
  typename uintd<D>::Type _dims = vector_to_uintd<D>( *(in->get_dimensions().get()) );
  typename intd<D>::Type dims;
  for( unsigned int i=0; i<D; i++ ){
    dims.vec[i] = (int)_dims.vec[i];
  }  
  
  if( D>3 ){
    std::cerr << std::endl << "partialDerivativeOperator::compute_laplace : internal error (only D<4 supported for now)." << std::endl;
    return -1;
  }

  _set_device();

  dim3 dimBlock( BLOCK_SIZE );
  dim3 dimGrid( prod(dims)/BLOCK_SIZE );
  
  // Invoke kernel

  mult_MH_kernel<REAL,T,D><<< dimGrid, dimBlock >>> (dims, in->get_data_ptr(), out->get_data_ptr(), _sigma->get_data_ptr() ,_norm->get_data_ptr());
  
  CHECK_FOR_CUDA_ERROR();

  _restore_device();
  
  return 0;
}

template< class REAL, class T, unsigned int D> int
cuVariableGaussOperator<REAL,T,D>::mult_MH_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate)
{
  cuNDArray<T> tmp;
  tmp.create(out->get_dimensions().get());

  int res = mult_M(in,&tmp);

  if (res < 0){
    return -1;
  }

  return mult_MH(&tmp,out,accumulate);
}

// Instantiations


template class EXPORTSOLVERS cuVariableGaussOperator<float, float, 1>;
template class EXPORTSOLVERS cuVariableGaussOperator<float, float, 2>;
template class EXPORTSOLVERS cuVariableGaussOperator<float, float, 3>;

/*
template class EXPORTSOLVERS cuVariableGaussOperator<float, float_complext::Type, 1>;
template class EXPORTSOLVERS cuVariableGaussOperator<float, float_complext::Type, 2>;
template class EXPORTSOLVERS cuVariableGaussOperator<float, float_complext::Type, 3>;
*/

template class EXPORTSOLVERS cuVariableGaussOperator<double, double, 1>;
template class EXPORTSOLVERS cuVariableGaussOperator<double, double, 2>;
template class EXPORTSOLVERS cuVariableGaussOperator<double, double, 3>;

/*
template class EXPORTSOLVERS cuVariableGaussOperator<double, double_complext::Type, 1>;
template class EXPORTSOLVERS cuVariableGaussOperator<double, double_complext::Type, 2>;
template class EXPORTSOLVERS cuVariableGaussOperator<double, double_complext::Type, 3>;
*/
