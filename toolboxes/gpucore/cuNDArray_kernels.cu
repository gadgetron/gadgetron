#include "cuNDArray.h"
#include "vector_td.h"

template <class T> 
__global__ void cuNDArray_permute_kernel(T* in, T* out, 
					 unsigned int ndim,
					 unsigned int* dims,
					 unsigned int* strides_out,
					 unsigned long int elements,
					 int shift_mode)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  unsigned long idx_out = 0;

  unsigned long idx_in_tmp = idx_in;
  if (idx_in < elements) {

    unsigned int cur_index;
    for (unsigned int i = 0; i < ndim; i++) {
      unsigned long idx_in_remainder = idx_in_tmp / dims[i];
      cur_index = idx_in_tmp-(idx_in_remainder*dims[i]); //cur_index = idx_in_tmp%dims[i];
      if (shift_mode < 0) { //IFFTSHIFT
	idx_out += ((cur_index+(dims[i]>>1))%dims[i])*strides_out[i];
      } else if (shift_mode > 0) { //FFTSHIFT
	idx_out += ((cur_index+((dims[i]+1)>>1))%dims[i])*strides_out[i];
      } else {
	idx_out += cur_index*strides_out[i];
      }
      idx_in_tmp = idx_in_remainder;
    }

    out[idx_in] = in[idx_out];

  }

}


template <class T> int cuNDArray_permute(cuNDArray<T>* in,
 					 cuNDArray<T>* out,
					 std::vector<unsigned int> order,
					 int shift_mode)
{
  cudaError_t err;

  T* in_ptr = in->data_;
  T* out_ptr = 0;

  if (out) {
    out_ptr = out->data_;
  } else {
    if (cudaMalloc((void**) &out_ptr, in->elements_*sizeof(T)) != cudaSuccess) {
      std::cerr << "cuNDArray_permute : Error allocating CUDA memory" << std::endl;
      out_ptr = 0;
      return -1;
    }
  }

  unsigned int* dims        = new unsigned int[in->dimensions_.size()];
  unsigned int* strides_out = new unsigned int[in->dimensions_.size()];
  if (!dims || !strides_out) {
    std::cerr << "cuNDArray_permute: failed to allocate temporary storage for arrays" << std::endl;
    return -1;
  }

  for (unsigned int i = 0; i < in->dimensions_.size(); i++) {
    dims[i] = in->dimensions_[order[i]];
    strides_out[i] = 1;
    
    for (unsigned int j = 0; j < order[i]; j++) {
      strides_out[i] *= in->dimensions_[j];
    }
  }

  unsigned int* dims_dev        = 0;
  unsigned int* strides_out_dev = 0;
  
  if (cudaMalloc((void**) &dims_dev, in->dimensions_.size()*sizeof(unsigned int)) != cudaSuccess) {
    std::cerr << "cuNDArray_permute : Error allocating CUDA dims memory" << std::endl;
    return -1;
  }
  
  if (cudaMalloc((void**) &strides_out_dev, in->dimensions_.size()*sizeof(unsigned int)) != cudaSuccess) {
    std::cerr << "cuNDArray_permute : Error allocating CUDA strides_out memory" << std::endl;
    return -1;
  }
  
  if (cudaMemcpy(dims_dev, dims, in->dimensions_.size()*sizeof(unsigned int), cudaMemcpyHostToDevice) !=
      cudaSuccess) {

    err = cudaGetLastError();
    std::cerr << "cuNDArray_permute : Error uploading dimensions to device, " 
	      << cudaGetErrorString(err) << std::endl;
    return -1;
  }

  if (cudaMemcpy(strides_out_dev, strides_out, in->dimensions_.size()*sizeof(unsigned int), cudaMemcpyHostToDevice) !=
      cudaSuccess) {
    std::cerr << "cuNDArray_permute : Error uploading strides to device" << std::endl;
    return -1;
  }

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)in->elements_/blockDim.x), 1, 1 );

  cuNDArray_permute_kernel<<< gridDim, blockDim >>>( in_ptr, out_ptr, in->dimensions_.size(), 
						     dims_dev, strides_out_dev, in->elements_, shift_mode);

  err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cuNDArray_permute : Error during kernel call: " << cudaGetErrorString(err) << std::endl;
    return -1;
  }

  if (cudaFree(dims_dev) != cudaSuccess) {
    err = cudaGetLastError();
    std::cerr << "cuNDArray_permute: failed to delete device memory (dims_dev) " 
	      << cudaGetErrorString(err) << std::endl;
    return -1;
  }

  if (cudaFree(strides_out_dev) != cudaSuccess) {
    err = cudaGetLastError();
    std::cerr << "cuNDArray_permute: failed to delete device memory (strides_out_dev) " 
	      << cudaGetErrorString(err) << std::endl;
    return -1;
  }
  
  delete [] dims;
  delete [] strides_out;

  if (!out) {
    std::vector<unsigned int> new_dims;
    for (unsigned int i = 0; i < in->dimensions_.size(); i++) {
      new_dims.push_back(in->dimensions_[order[i]]);
    }
    in->dimensions_ = new_dims;
    if (cudaFree(in->data_) != cudaSuccess) {
	std::cerr << "cuNDArray_permute: failed to delete device memory" << std::endl;
	return -1;
    }
    in->data_ = out_ptr;
  }

  return 0;
}

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<int>* in,
				 cuNDArray<int>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<int2>* in,
				 cuNDArray<int2>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<int3>* in,
				 cuNDArray<int3>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<int4>* in,
				 cuNDArray<int4>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<unsigned int>* in,
				 cuNDArray<unsigned int>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uint2>* in,
				 cuNDArray<uint2>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uint3>* in,
				 cuNDArray<uint3>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uint4>* in,
				 cuNDArray<uint4>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<float>* in,
				 cuNDArray<float>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<float2>* in,
				 cuNDArray<float2>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<float3>* in,
				 cuNDArray<float3>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<float4>* in,
				 cuNDArray<float4>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<double>* in,
				 cuNDArray<double>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<double2>* in,
				 cuNDArray<double2>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<double3>* in,
				 cuNDArray<double3>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<double4>* in,
				 cuNDArray<double4>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<intd<1>::Type>* in,
				 cuNDArray<intd<1>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<intd<2>::Type>* in,
				 cuNDArray<intd<2>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<intd<3>::Type>* in,
				 cuNDArray<intd<3>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<intd<4>::Type>* in,
				 cuNDArray<intd<4>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uintd<1>::Type>* in,
				 cuNDArray<uintd<1>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uintd<2>::Type>* in,
				 cuNDArray<uintd<2>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uintd<3>::Type>* in,
				 cuNDArray<uintd<3>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uintd<4>::Type>* in,
				 cuNDArray<uintd<4>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<floatd<1>::Type>* in,
				 cuNDArray<floatd<1>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<floatd<2>::Type>* in,
				 cuNDArray<floatd<2>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<floatd<3>::Type>* in,
				 cuNDArray<floatd<3>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<floatd<4>::Type>* in,
				 cuNDArray<floatd<4>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<doubled<1>::Type>* in,
				 cuNDArray<doubled<1>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<doubled<2>::Type>* in,
				 cuNDArray<doubled<2>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<doubled<3>::Type>* in,
				 cuNDArray<doubled<3>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<doubled<4>::Type>* in,
				 cuNDArray<doubled<4>::Type>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<intd1>* in,
				 cuNDArray<intd1>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<intd2>* in,
				 cuNDArray<intd2>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<intd3>* in,
				 cuNDArray<intd3>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<intd4>* in,
				 cuNDArray<intd4>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uintd1>* in,
				 cuNDArray<uintd1>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uintd2>* in,
				 cuNDArray<uintd2>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uintd3>* in,
				 cuNDArray<uintd3>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<uintd4>* in,
				 cuNDArray<uintd4>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<floatd1>* in,
				 cuNDArray<floatd1>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<floatd2>* in,
				 cuNDArray<floatd2>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<floatd3>* in,
				 cuNDArray<floatd3>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<floatd4>* in,
				 cuNDArray<floatd4>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<doubled1>* in,
				 cuNDArray<doubled1>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<doubled2>* in,
				 cuNDArray<doubled2>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<doubled3>* in,
				 cuNDArray<doubled3>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<doubled4>* in,
				 cuNDArray<doubled4>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);
				   
template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<float_complext>* in,
				 cuNDArray<float_complext>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);

template EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<double_complext>* in,
				 cuNDArray<double_complext>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);
