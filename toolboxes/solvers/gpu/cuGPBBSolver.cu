#include "cuGPBBSolver.h"
#include "complext.h"

#define MAX_THREADS_PER_BLOCK 512

using namespace Gadgetron;
template <class T> __global__ void filter_kernel(T* x, T* g, int elements){
  typedef typename realType<T>::Type REAL;
  const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < elements){
    if ( real(x[idx]) <= REAL(0) && real(g[idx]) > 0) g[idx]=T(0);
  }
}

template <class T> void Gadgetron::cuGPBBSolver<T>::solver_non_negativity_filter(Gadgetron::cuNDArray<T>* x , Gadgetron::cuNDArray<T>* g){
  int elements = g->get_number_of_elements();

  int threadsPerBlock = std::min(elements,MAX_THREADS_PER_BLOCK);
  dim3 dimBlock( threadsPerBlock);
  int totalBlocksPerGrid = std::max(1,elements/MAX_THREADS_PER_BLOCK);
  dim3 dimGrid(totalBlocksPerGrid);

  filter_kernel<T><<<dimGrid,dimBlock>>>(x->get_data_ptr(),g->get_data_ptr(),elements);
}

template<typename T>
struct negate2 : public thrust::unary_function<T,T>
{
  __host__ __device__ T operator()(const T &x) const {return -x;}
};

template<typename REAL, typename T>
struct reciprocal_clamp_functor : public thrust::unary_function<T,T>
{
  const REAL clamp;
  reciprocal_clamp_functor(REAL _clamp) : clamp(_clamp){}
  __host__ __device__
  T operator()(const T &x) const{
    if (real(x) < clamp) return T(0);
    else return T(1)/x;
  }
};

template <class T> void Gadgetron::cuGPBBSolver<T>::solver_reciprocal_clamp( Gadgetron::cuNDArray<T>* x,REAL threshold) 
{
  thrust::device_ptr<T> dev_ptr(x->get_data_ptr());
  thrust::transform(dev_ptr, dev_ptr + x->get_number_of_elements(),dev_ptr,reciprocal_clamp_functor<REAL,T>(threshold));
};

template class EXPORTGPUSOLVERS Gadgetron::cuGPBBSolver<float>;
template class EXPORTGPUSOLVERS Gadgetron::cuGPBBSolver<double>;
