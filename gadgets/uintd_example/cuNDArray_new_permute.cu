#include "cuNDArray_new_permute.h"

#include <cuComplex.h>

#include "uintd.h"

template <int ASIZE>  __host__ __device__ void nidx_to_co(unsigned int& idx, uintd<ASIZE>& dims, uintd<ASIZE>& co)
{
  unsigned int idx_tmp = idx;
  for (unsigned int i = 0; i < ASIZE; i++) {
    co.d[i] = idx_tmp%dims.d[i];
    idx_tmp -= co.d[i];
    idx_tmp /= dims.d[i];
  }
} 

template <int ASIZE>  __host__ __device__ unsigned int nco_to_idx(uintd<ASIZE>& dims, uintd<ASIZE>& co)
{
  unsigned int idx = 0;
  unsigned long block_size = 1;
  for (unsigned int i = 0; i < ASIZE; i++) {
    idx += block_size*co.d[i];
     block_size *= dims.d[i];
  }
  return idx;
} 

template <int ASIZE>  __host__ __device__ unsigned int nco_to_idx(uintd<ASIZE>& dims, uintd<ASIZE>& co, uintd<ASIZE>& order)
{
  unsigned int idx = 0;
  unsigned long block_size = 1;
  for (unsigned int i = 0; i < ASIZE; i++) {
    idx += block_size*co.d[ order.d[i] ];
     block_size *= dims.d[ order.d[i] ];
  }
  return idx;
} 


template <class T, int DIMENSIONS> __global__ void cuNDArray_permute_kernel_new(T* in, T* out, 
								       uintd<DIMENSIONS> dims,
								       uintd<DIMENSIONS> order,
								       unsigned long int elements,
								       int shift_mode)
{
  unsigned int idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  unsigned int idx_out = 0;

  if (idx_in < elements) {
    uintd<DIMENSIONS> co;
    nidx_to_co(idx_in, dims, co);
    idx_out = nco_to_idx(co, dims, order);
    //out[idx_out] = in[idx_in];
    out[idx_in].x = idx_out;//in[idx_in];
  }
}


template < class T, int DIMENSIONS > int cuNDArray_new_permute_dims(cuNDArray<T>* in, cuNDArray<T>* out, std::vector<unsigned int> order)
{
  cudaError_t err;

  T* in_ptr = in->get_data_ptr();
  T* out_ptr = 0;

  if (out) {
    out_ptr = out->get_data_ptr();
  } else {
    if (cudaMalloc((void**) &out_ptr, in->get_number_of_elements()*sizeof(T)) != cudaSuccess) {
      std::cerr << "cuNDArray_permute : Error allocating CUDA memory" << std::endl;
      out_ptr = 0;
      return -1;
    }
  }

  std::vector<unsigned int> dimensions = in->get_dimensions();
  uintd<DIMENSIONS> dims_tmp;
  uintd<DIMENSIONS> order_tmp;

  for (unsigned int i = 0; i < dimensions.size(); i++) {
    dims_tmp[i] = dimensions[i];
    order_tmp[i] = order[i];
  }

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)in->get_number_of_elements()/blockDim.x), 1, 1 );
  cuNDArray_permute_kernel_new<<< gridDim, blockDim >>>( in_ptr, out_ptr, dims_tmp, order_tmp, in->get_number_of_elements(), 0);

  err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cuNDArray_permute : Error during kernel call: " << cudaGetErrorString(err) << std::endl;
    return -1;
  }

  if (!out) {
    std::vector<unsigned int> new_dims;
    for (unsigned int i = 0; i < dimensions.size(); i++) {
      new_dims.push_back(dimensions[order[i]]);
    }
    in->create(new_dims, out_ptr, true);
  }

  return 0;
}


template <class T> int cuNDArray_new_permute(cuNDArray<T>* in, cuNDArray<T>* out, std::vector<unsigned int> order)
{
  unsigned int ndim = in->get_number_of_dimensions();

  switch (ndim) {
  case 1:
    return cuNDArray_new_permute_dims<T, 1>(in, out, order); 
  case 2:
    return cuNDArray_new_permute_dims<T, 2>(in, out, order); 
  case 3:
    return cuNDArray_new_permute_dims<T, 3>(in, out, order); 
  case 4:
    return cuNDArray_new_permute_dims<T, 4>(in, out, order); 
  default:
    std::cerr << "Permute not supported for number of dimensions = " << ndim << std::endl;
    return -1;
  }
}

template int cuNDArray_new_permute(cuNDArray<float2>* in, cuNDArray<float2>* out, std::vector<unsigned int> order);
//template int cuNDArray_new_permute(cuNDArray<float>* in, cuNDArray<float>* out, std::vector<unsigned int> order);

/*


template void nidx_to_co(unsigned int& idx, uintd<1>& dims, uintd<1>& co);
template void nidx_to_co(unsigned int& idx, uintd<2>& dims, uintd<2>& co);
template void nidx_to_co(unsigned int& idx, uintd<3>& dims, uintd<3>& co);
template void nidx_to_co(unsigned int& idx, uintd<4>& dims, uintd<4>& co);

template unsigned int nco_to_idx(uintd<1>& dims, uintd<1>& co);
template unsigned int nco_to_idx(uintd<2>& dims, uintd<2>& co);
template unsigned int nco_to_idx(uintd<3>& dims, uintd<3>& co);
template unsigned int nco_to_idx(uintd<4>& dims, uintd<4>& co);

template unsigned int nco_to_idx(uintd<1>& dims, uintd<1>& co, uintd<1>& order);
template unsigned int nco_to_idx(uintd<2>& dims, uintd<2>& co, uintd<2>& order);
template unsigned int nco_to_idx(uintd<3>& dims, uintd<3>& co, uintd<3>& order);
template unsigned int nco_to_idx(uintd<4>& dims, uintd<4>& co, uintd<4>& order);
*/
