#include "cuGpBbSolver.h"
#include "complext.h"

#define MAX_THREADS_PER_BLOCK 512

using namespace Gadgetron;
template <class T> __global__ void filter_kernel(T* x, T* g, int elements){
  const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < elements){
    if ( x[idx] <= T(0) && g[idx] > 0) g[idx]=T(0);
  }
}

template <class REAL> __global__ void filter_kernel(complext<REAL>* x, complext<REAL>* g, int elements){
  const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < elements){
    if ( real(x[idx]) <= REAL(0) && real(g[idx]) > 0) g[idx].vec[0] = REAL(0);
    g[idx].vec[1]=REAL(0);
  }
}

template <class T> void Gadgetron::cuGpBbSolver<T>::
solver_non_negativity_filter(Gadgetron::cuNDArray<T>* x , Gadgetron::cuNDArray<T>* g)
{
  int elements = g->get_number_of_elements();

  int threadsPerBlock = std::min(elements,MAX_THREADS_PER_BLOCK);
  dim3 dimBlock( threadsPerBlock);
  int totalBlocksPerGrid = std::max(1,elements/MAX_THREADS_PER_BLOCK);
  dim3 dimGrid(totalBlocksPerGrid);

  filter_kernel<typename realType<T>::Type><<<dimGrid,dimBlock>>>(x->get_data_ptr(),g->get_data_ptr(),elements);
}


template class EXPORTGPUSOLVERS Gadgetron::cuGpBbSolver<float>;
template class EXPORTGPUSOLVERS Gadgetron::cuGpBbSolver<double>;

template class EXPORTGPUSOLVERS Gadgetron::cuGpBbSolver< complext<float> >;
template class EXPORTGPUSOLVERS Gadgetron::cuGpBbSolver< complext<double> >;
