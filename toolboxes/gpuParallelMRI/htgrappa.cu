#include "htgrappa.h"
#include <cublas.h>

#include <cublas_v2.h>
#include "cposv_wrapper.h"
#include "hoNDArray_fileio.h"
#include "cuNDFFT.h"
#include "GPUTimer.h"

int2 vec_to_int2(std::vector<unsigned int> vec)
{
  int2 ret; ret.x = 0; ret.y = 0;
  if (vec.size() < 2) {
    std::cout << "vec_to_uint2 dimensions of vector too small" << std::endl;
    return ret;
  }

  ret.x = vec[0]; ret.y = vec[1];
  return ret;
}


__global__ void clear_array(float2* in, unsigned long int elements)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < elements) {
    in[idx_in].x = 0.0;
    in[idx_in].y = 0.0;
  }
}

int clear(cuNDArray<float2>* in)
{
  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)in->get_number_of_elements()/blockDim.x), 1, 1 );

  clear_array<<< gridDim, blockDim >>>( in->get_data_ptr(), in->get_number_of_elements());

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "clear : Error during kernel call: " << cudaGetErrorString(err) << std::endl;
    return -1;
  }

  return 0;
}

template <class T> int write_cuNDArray_to_disk(cuNDArray<T>* a, const char* filename)
{
  hoNDArray<T> host = a->to_host();
  write_nd_array<cuFloatComplex>(host, filename);
  return 0;
}

template <class T> __global__ void form_grappa_system_matrix_kernel_2d(T* ref_data,
								       int2 dims,
								       int coils,
								       int2 ros,
								       int2 ros_offset,
								       int2 kernel_size,
								       int acceleration_factor,
								       int set_number,			       
								       T* out_matrix,
								       T* b)
{
   long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
   int klocations = ros.x*ros.y;
   int image_elements = dims.x*dims.y;
   //int coefficients = kernel_size.x*kernel_size.y*coils;
   if (idx_in < klocations) {
     //unsigned int y = idx_in/ros.x;
     //unsigned int x = idx_in - y*ros.x;
     unsigned int x = idx_in/ros.y;
     unsigned int y = idx_in - x*ros.y;
     unsigned int idx_ref = 0;
     unsigned int coeff_counter = 0;

     int kernel_size_x = kernel_size.x;
     int kernel_size_y = kernel_size.y;
     for (int c = 0; c < coils; c++) {
	 for (int ky = -((kernel_size_y*acceleration_factor)>>1)+set_number+1; 
	      ky < ((kernel_size_y*acceleration_factor+1)>>1); ky+=acceleration_factor) {
	   for (int kx = -(kernel_size_x>>1); kx < ((kernel_size_x+1)>>1); kx++) {
	     idx_ref = c*image_elements + x+kx+ros_offset.x + (y+ky+ros_offset.y)*dims.x;
	     //out_matrix[idx_in*coefficients+coeff_counter++] = ref_data[idx_ref];
	     out_matrix[idx_in+(coeff_counter++)*klocations] = ref_data[idx_ref];
	   
	 }
       }
     }

     for (unsigned int c = 0; c < coils; c++) {
       //b[idx_in*coils + c] = ref_data[c*image_elements + y*dims.x+x];
       b[idx_in + c*klocations] = ref_data[c*image_elements + (y+ros_offset.y)*dims.x+(x+ros_offset.x)];
     }
   }
}

template <class T> __global__ void copy_grappa_coefficients_to_kernel_2d(T* coeffs,
									 T* kernel,
									 int coils,
									 int2 kernel_size,
									 int acceleration_factor,
									 int set)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  
  unsigned int coefficients_in_set = coils*kernel_size.x*kernel_size.y*coils;

  if (idx_in < coefficients_in_set) {
    int idx_in_tmp = idx_in;
    int kx = idx_in%kernel_size.x;
    idx_in = (idx_in-kx)/kernel_size.x;
    int ky = idx_in%kernel_size.y;
    idx_in = (idx_in-ky)/kernel_size.y;
    int coil = idx_in%coils;
    idx_in = (idx_in-coil)/coils;
    int coilg = idx_in;

    kernel[coilg*coils*(kernel_size.y*acceleration_factor)*kernel_size.x +
	   coil*(kernel_size.y*acceleration_factor)*kernel_size.x +
	   (ky*acceleration_factor + set + 1)*kernel_size.x + kx] = coeffs[idx_in_tmp];

    if ((coil == coilg) && (kx == 0) && (ky == 0) && (set == 0)) {
      kernel[coilg*coils*(kernel_size.y*acceleration_factor)*kernel_size.x +
	     coil*(kernel_size.y*acceleration_factor)*kernel_size.x +
	     ((kernel_size.y>>1)*acceleration_factor)*kernel_size.x + (kernel_size.x>>1) ].x = 1;
      
    }
  }
}

template <class T> __global__ void copy_grappa_kernel_to_kspace_2d(T* kernel,
								   T* out,
								   int2 dims,
								   int2 kernel_size,
								   int coils)

{
  
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;

  if (idx_in < kernel_size.x*kernel_size.y*coils) {
    int idx_in_tmp = idx_in;
    int kx = idx_in%kernel_size.x;
    idx_in = (idx_in-kx)/kernel_size.x;
    int ky = idx_in%kernel_size.y;
    idx_in = (idx_in-ky)/kernel_size.y;
    int coil = idx_in;

    int outx = -(kx- (kernel_size.x>>1)) + (dims.x>>1); //Flipping the kernel for conv
    int outy = -(ky- (kernel_size.y>>1)) + (dims.y>>1);

    out[coil*dims.x*dims.y + outy*dims.x + outx] = kernel[idx_in_tmp];
  }

}

__global__ void scale_and_add_unmixing_coeffs(cuFloatComplex* unmixing,
								 cuFloatComplex* csm,
								 cuFloatComplex* out,
								 int elements,
								 int coils)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;

  cuFloatComplex tmp;
  if (idx_in < elements) {
    for (int c = 0; c < coils; c++) {
      tmp = cuCmulf(unmixing[c*elements + idx_in],cuConjf(csm[idx_in])); 
      out[c*elements + idx_in].x += tmp.x;
      out[c*elements + idx_in].y += tmp.y;
    }
  }
}

template <class T> int htgrappa_calculate_grappa_unmixing(cuNDArray<T>* ref_data, 
							  cuNDArray<T>* b1,
							  unsigned int acceleration_factor,
							  std::vector<unsigned int> kernel_size,
							  cuNDArray<T>* out_mixing_coeff)
{

  if (!ref_data->dimensions_equal(*b1) ||
      !ref_data->dimensions_equal(*out_mixing_coeff)) {
    std::cerr << "htgrappa_calculate_grappa_unmixing: Dimensions mismatch" << std::endl;
    return -1;
  }

  if (kernel_size.size() != (ref_data->get_number_of_dimensions()-1)) {
    std::cerr << "htgrappa_calculate_grappa_unmixing: Kernel size does not match the data dimensions" << std::endl;
    return -1;
  }

  if (ref_data->get_number_of_dimensions() > 3) {
    std::cerr << "htgrappa_calculate_grappa_unmixing: Not yet implemented for 3D" << std::endl;
    return -1;
  }
    
  //Calculate region of support + offsets
  std::vector<unsigned int> ros = ref_data->get_dimensions();
  ros.pop_back(); //Remove the number of coils

  std::vector<unsigned int> ros_offset(ref_data->get_number_of_dimensions(),0);
  unsigned long int kspace_locations = 1;
  for (unsigned int i = 0; i < ros.size(); i++) {
    if (i > 0) {
      ros[i] -= (kernel_size[i]*acceleration_factor);
    } else {
      ros[i] -= kernel_size[i];
    }
    ros_offset[i] = (ref_data->get_size(i)-ros[i])>>1;
    kspace_locations *= ros[i];
  } 

  unsigned int coils = ref_data->get_size(ref_data->get_number_of_dimensions()-1);

  std::vector<unsigned int> sys_matrix_size; 
  sys_matrix_size.push_back(kspace_locations);
  sys_matrix_size.push_back(coils*kernel_size[0]*kernel_size[1]);
  
  std::vector<unsigned int> b_size;
  b_size.push_back(kspace_locations);
  b_size.push_back(coils);

  cuNDArray<T> system_matrix;
  if (!system_matrix.create(sys_matrix_size)) {
    std::cout << "htgrappa_calculate_grappa_unmixing: Unable to allocate device memory for system matrix" << std::endl;
    return -1;
  }

  clear(&system_matrix);

  cuNDArray<T> b;
  if (!b.create(b_size)) {
    std::cout << "htgrappa_calculate_grappa_unmixing: Unable to allocate device memory for right hand sides" << std::endl;
    return -1;
  }

  int2 dims = vec_to_int2(ref_data->get_dimensions());
  int2 dros = vec_to_int2(ros);
  int2 dros_offset = vec_to_int2(ros_offset);
  int2 dkernel_size = vec_to_int2(kernel_size);

  int n = coils*kernel_size[0]*kernel_size[1];
  int m = kspace_locations;

  std::vector<unsigned int> AHA_dims(2,n);
  cuNDArray<T> AHA;
  if (!AHA.create(AHA_dims)) {
    std::cout << "htgrappa_calculate_grappa_unmixing: Unable to allocate device memory for AHA" << std::endl;
    return -1;
  }

  std::vector<unsigned int> AHrhs_dims;
  AHrhs_dims.push_back(n);
  AHrhs_dims.push_back(coils);

  cuNDArray<T> AHrhs;
  if (!AHrhs.create(AHrhs_dims)) {
    std::cout << "htgrappa_calculate_grappa_unmixing: Unable to allocate device memory for AHrhs" << std::endl;
    return -1;
  }


  cublasHandle_t handle;
  if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS) {
    std::cerr << "htgrappa_calculate_grappa_unmixing: unable to create cublas handle" << std::endl;
    return -1;

  }

  std::vector<unsigned int> gkernel_dims;
  gkernel_dims.push_back(kernel_size[0]);
  gkernel_dims.push_back(kernel_size[1]*acceleration_factor);
  gkernel_dims.push_back(coils);
  gkernel_dims.push_back(coils);
  cuNDArray<T> gkernel;
  if (!gkernel.create(gkernel_dims)) {
    std::cerr << "htgrappa_calculate_grappa_unmixing: Unable to allocate array for GRAPPA kernel" << std::endl;
    return -1;
  }

  clear(&gkernel);
  
  for (unsigned int set = 0; set < acceleration_factor-1; set++) {
    //std::cout << "Calculating coefficients for set " << set << std::endl;

    //std::cout << "dros.x = " << dros.x << ", dros.y = " << dros.y << std::endl;

    dim3 blockDim(512,1,1);
    dim3 gridDim((unsigned int) ceil((1.0f*kspace_locations)/blockDim.x), 1, 1 );
  
    form_grappa_system_matrix_kernel_2d<<< gridDim, blockDim >>>( ref_data->get_data_ptr(), dims, coils, dros, dros_offset,
								  dkernel_size, acceleration_factor, set, 
								  system_matrix.get_data_ptr(),
								  b.get_data_ptr());

    cudaError_t err = cudaGetLastError();
    if( err != cudaSuccess ){
      std::cerr << "htgrappa_calculate_grappa_unmixing: Unable to form system matrix: " << 
	cudaGetErrorString(err) << std::endl;
      return -1;
    }

    //write_cuNDArray_to_disk(&system_matrix,"A.cplx");
    //write_cuNDArray_to_disk(&b,"b.cplx");

    cuFloatComplex alpha = make_float2(1.0,0.0);
    cuFloatComplex beta = make_float2(0.0, 0.0);

    cublasStatus_t stat = cublasCgemm(handle, CUBLAS_OP_C, CUBLAS_OP_N,
				      n,n,m,&alpha,
				      system_matrix.get_data_ptr(), m,
				      system_matrix.get_data_ptr(), m,
				      &beta, AHA.get_data_ptr(), n);
    
    if (stat != CUBLAS_STATUS_SUCCESS) {
      std::cerr << "htgrappa_calculate_grappa_unmixing: Failed to form AHA product using cublas gemm" << std::endl;
      std::cerr << "---- cublas error code " << stat << std::endl;
      return -1;
    }

    //write_cuNDArray_to_disk(&AHA,"AHA.cplx");

    {
      //GPUTimer timer("GRAPPA cublas gemm");
      stat = cublasCgemm(handle, CUBLAS_OP_C, CUBLAS_OP_N,
		       n,coils,m,&alpha,
		       system_matrix.get_data_ptr(), m,
		       b.get_data_ptr(), m,
		       &beta, AHrhs.get_data_ptr(), n);
    
    }
    //write_cuNDArray_to_disk(&AHrhs,"AHrhs.cplx");
 
    if (stat != CUBLAS_STATUS_SUCCESS) {
      std::cerr << "htgrappa_calculate_grappa_unmixing: Failed to form AHA product using cublas gemm" << std::endl;
      std::cerr << "---- cublas error code " << stat << std::endl;
      return -1;
    }

    if (cposv_wrapper(&AHA, &AHrhs) < 0) {
      std::cerr << "htgrappa_calculate_grappa_unmixing: Error calling cgels" << std::endl;
      return -1;
    }
    
    //write_cuNDArray_to_disk(&AHrhs,"AHrhs_solution.cplx");

    gridDim = dim3((unsigned int) ceil((1.0f*n*coils)/blockDim.x), 1, 1 );
  
    copy_grappa_coefficients_to_kernel_2d<<< gridDim, blockDim >>>( AHrhs.get_data_ptr(), 
								    gkernel.get_data_ptr(),
								    coils,
								    dkernel_size,
								    acceleration_factor,
								    set);

    //write_cuNDArray_to_disk(&gkernel,"kernel.cplx");
 
    err = cudaGetLastError();
    if( err != cudaSuccess ){
      std::cerr << "htgrappa_calculate_grappa_unmixing: Failed to copy calculated coefficients to kernel: " << 
	cudaGetErrorString(err) << std::endl;
      return -1;
    }

  }

  cuNDArray<T> tmp_mixing;
  if (!tmp_mixing.create(out_mixing_coeff->get_dimensions())) {
    std::cerr << "htgrappa_calculate_grappa_unmixing: Unable to create temp mixing storage on device." << std::endl;
    return -1;
  }
    

  int kernel_elements = gkernel.get_number_of_elements()/coils;
  int total_elements = tmp_mixing.get_number_of_elements()/coils;
  dkernel_size.y *= acceleration_factor;
  cuNDFFT ft;
  std::vector<unsigned int> ft_dims(2,0);ft_dims[1] = 1;
  clear(out_mixing_coeff);
  for (unsigned int c = 0; c < coils; c++) {
    clear(&tmp_mixing);

    dim3 blockDim(512,1,1);
    dim3 gridDim((unsigned int) ceil((1.0f*kernel_elements)/blockDim.x), 1, 1 ); 
    copy_grappa_kernel_to_kspace_2d<<< gridDim, blockDim >>>((gkernel.get_data_ptr()+(c*kernel_elements)),
							     tmp_mixing.get_data_ptr(),
							     dims,
							     dkernel_size,
							     coils);

    cudaError_t err = cudaGetLastError();
    if( err != cudaSuccess ){
      std::cerr << "htgrappa_calculate_grappa_unmixing: Unable to pad GRAPPA kernel: " <<
	cudaGetErrorString(err) << std::endl;
      return -1;
    }

    ft.ifft(&tmp_mixing,ft_dims);

    gridDim = dim3((unsigned int) ceil(1.0f*total_elements/blockDim.x), 1, 1 ); 
    scale_and_add_unmixing_coeffs<<< gridDim, blockDim >>>(tmp_mixing.get_data_ptr(),
							   (b1->get_data_ptr()+ c*total_elements),
							   out_mixing_coeff->get_data_ptr(),
							   total_elements,
							   coils);
    err = cudaGetLastError();
    if( err != cudaSuccess ){
      std::cerr << "htgrappa_calculate_grappa_unmixing: scale and add mixing coeffs: " << 
	cudaGetErrorString(err) << std::endl;
      return -1;
    }
  }

  cublasDestroy(handle);
 
  return 0;
}



//Template instanciation
template int htgrappa_calculate_grappa_unmixing(cuNDArray<float2>* ref_data, 
						cuNDArray<float2>* b1,
						unsigned int acceleration_factor,
						std::vector<unsigned int> kernel_size,
						cuNDArray<float2>* out_mixing_coeff);
