#include "cuTvOperator.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "vector_td_utilities.h"
#include "complext.h"
#include <iostream>
#include "check_CUDA.h"
#include "cudaDeviceManager.h"

using namespace Gadgetron;

template<class REAL, class T, unsigned int D> static inline  __device__ REAL gradient(const T* in, const vector_td<int,D>& dims, vector_td<int,D>& co)
{
  REAL grad = REAL(0);
  T xi = in[co_to_idx<D>((co+dims)%dims,dims)];
  for (int i = 0; i < D; i++){
    co[i]+=1;
    T dt = in[co_to_idx<D>((co+dims)%dims,dims)];
    grad += norm(xi-dt);
    co[i]-=1;
  }
  return sqrt(grad);
}


template<class REAL, class T, unsigned int D> static __global__ void tvGradient_kernel(const T* in, T* out, const vector_td<int,D> dims,REAL limit,REAL weight)
{
  const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < prod(dims) ){
    T xi = in[idx];
    T result=T(0);

    vector_td<int,D> co = idx_to_co<D>(idx, dims);

    REAL grad = gradient<REAL,T,D>(in,dims,co);

    if (grad > limit) {
      result += REAL(D)*xi/grad;
      for (int i = 0; i < D; i++){
	co[i]+=1;
	result -= in[co_to_idx<D>((co+dims)%dims,dims)]/grad;
	co[i]-=1;
      }
    }

    for (int i = 0; i < D; i++){
      co[i]-=1;
      grad = gradient<REAL,T,D>(in,dims,co);
      if (grad > limit) {
	result +=(xi-in[co_to_idx<D>((co+dims)%dims,dims)])/grad;
      }
      co[i]+=1;
    }
    out[idx] += result*weight;
  }
}





template<class T, unsigned int D> void cuTvOperator<T,D>::gradient (cuNDArray<T> * in,cuNDArray<T> * out, bool accumulate)
{
  if (!accumulate) 
    clear(out);

  const typename intd<D>::Type dims = to_intd( from_std_vector<size_t,D>(*(in->get_dimensions())));
  int elements = in->get_number_of_elements();

  int threadsPerBlock =std::min(prod(dims),cudaDeviceManager::Instance()->max_blockdim());
  dim3 dimBlock( threadsPerBlock);
  int totalBlocksPerGrid = std::max(1,prod(dims)/cudaDeviceManager::Instance()->max_blockdim());
  dim3 dimGrid(totalBlocksPerGrid);

  for (int i =0; i < (elements/prod(dims)); i++){
    tvGradient_kernel<<<dimGrid,dimBlock>>>(in->get_data_ptr()+i*prod(dims),out->get_data_ptr()+i*prod(dims),dims,limit_,this->weight_);
  }
  
  cudaDeviceSynchronize();
  CHECK_FOR_CUDA_ERROR();
}

template<class REAL, class T, unsigned int D> static __global__ void tvMagnitude_kernel(const T* in,T* out,const vector_td<int,D> dims,REAL limit,REAL weight)
{
  const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < prod(dims) ){
    vector_td<int,D> co = idx_to_co<D>(idx, dims);
    REAL grad = gradient<REAL,T,D>(in,dims,co);
    out[idx] = grad*weight;
  }
}


template<class T, unsigned int D> typename realType<T>::Type cuTvOperator<T,D>::magnitude (cuNDArray<T> * in)
{
  cuNDArray<T> out(in->get_dimensions());
  const typename intd<D>::Type dims = to_intd( from_std_vector<size_t,D>(*(in->get_dimensions())));
  int elements = in->get_number_of_elements();

  int threadsPerBlock =std::min(prod(dims),cudaDeviceManager::Instance()->max_blockdim());
  dim3 dimBlock( threadsPerBlock);
  int totalBlocksPerGrid = std::max(1,prod(dims)/cudaDeviceManager::Instance()->max_blockdim());
  dim3 dimGrid(totalBlocksPerGrid);

  for (int i =0; i < (elements/prod(dims)); i++){
    tvMagnitude_kernel<<<dimGrid,dimBlock>>>(in->get_data_ptr()+i*prod(dims),out.get_data_ptr()+i*prod(dims),dims,limit_,this->weight_);
  }

  cudaDeviceSynchronize();
  CHECK_FOR_CUDA_ERROR();
  return asum(&out);
}

template class EXPORTGPUOPERATORS cuTvOperator<float,1>;
template class EXPORTGPUOPERATORS cuTvOperator<float,2>;
template class EXPORTGPUOPERATORS cuTvOperator<float,3>;
template class EXPORTGPUOPERATORS cuTvOperator<float,4>;

template class EXPORTGPUOPERATORS cuTvOperator<double,1>;
template class EXPORTGPUOPERATORS cuTvOperator<double,2>;
template class EXPORTGPUOPERATORS cuTvOperator<double,3>;
template class EXPORTGPUOPERATORS cuTvOperator<double,4>;

template class EXPORTGPUOPERATORS cuTvOperator<float_complext,1>;
template class EXPORTGPUOPERATORS cuTvOperator<float_complext,2>;
template class EXPORTGPUOPERATORS cuTvOperator<float_complext,3>;
template class EXPORTGPUOPERATORS cuTvOperator<float_complext,4>;

template class EXPORTGPUOPERATORS cuTvOperator<double_complext,1>;
template class EXPORTGPUOPERATORS cuTvOperator<double_complext,2>;
template class EXPORTGPUOPERATORS cuTvOperator<double_complext,3>;
template class EXPORTGPUOPERATORS cuTvOperator<double_complext,4>;
