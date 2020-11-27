#include "spirit_calibration.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_reductions.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_blas.h"
#include "cuNDFFT.h"
#include "cudaDeviceManager.h"
#include "setup_grid.h"
#include "complext.h"
#include "CUBLASContextProvider.h"
#include "GPUTimer.h"
#include "hoNDArray_fileio.h"
#include "htgrappa.h"

#include <cublas_v2.h>
//#include <cula_lapack_device.h>

namespace Gadgetron {

  static __global__ void 
  compute_system_matrix_kernel( intd2 dims,
                                int num_coils,
                                int kernel_size,
                                float_complext *kspace,
                                float_complext *A )
  {
    // The grid contains one thread per coil element. 
    // Each thread reads its corresponding data element and is responsible 
    // for filling into the corresponding kernel_size*kernel entries in the matrix.
    //
    // The storage format is column major due to BLAS/LAPACK conventions.
    // This increases the overhead of writes in this kernel (they are non-coaslesced and MANY). 
    // TODO: optimize for performance.
    //
    
    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const int elements_per_coil = prod(dims);

    if( idx < elements_per_coil*num_coils ){

      // Get the k-space value for this thread
      //

      float_complext val = kspace[idx];

      const int num_kernel_elements = kernel_size*kernel_size-1;
      const int coil = idx/elements_per_coil;
      const int idx_in_coil = idx-coil*elements_per_coil;

      // Loop over the number of outputs produced per thread
      //

      const int half_kernel_size = kernel_size>>1;

      for( int j = -half_kernel_size; j<half_kernel_size+1; j++ ){ // row iterator
        for( int i = -half_kernel_size; i<half_kernel_size+1; i++ ){ // column iterator

          if( j==0 && i==0 ) continue; // The weight of the central points is set to 0

          int kernel_idx = co_to_idx( intd2(i+half_kernel_size,j+half_kernel_size), intd2(kernel_size,kernel_size) );
          if( (j==0 && i>0) || j>0 ) kernel_idx--;

          const int m = 
            (idx_in_coil+j*dims[0]+i+elements_per_coil)%elements_per_coil; // row idx

          const int n = 
            coil*num_kernel_elements + kernel_idx;
          
          const int A_idx = 
            n*elements_per_coil + m; // Column major storage

          A[A_idx] = val;
        }
      }      
    }
  }
  

  static __global__ void 
  write_convolution_masks_kernel( intd2 dims,
                                  int num_coils,
                                  int kernel_size,
                                  float_complext *kernels,
                                  float_complext *kspace )
  {
    // Write out convolution masks in the center of kspace
    // - thus prepare for FFT into image space
    //

    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const int elements_per_coil = prod(dims);

    if( idx < elements_per_coil*num_coils*num_coils ){

      const int half_kernel_size = kernel_size>>1;
      const int num_kernel_elements = kernel_size*kernel_size-1;
      const int batch = idx/(elements_per_coil*num_coils);
      const int idx_in_batch = idx-batch*elements_per_coil*num_coils;
      const int coil = idx_in_batch/elements_per_coil;
      const int idx_in_coil = idx_in_batch-coil*elements_per_coil;
      const intd2 co = idx_to_co( idx_in_coil, dims ) - (dims>>1);

      if( co[1] >= -half_kernel_size && co[1] <= half_kernel_size && 
          co[0] >= -half_kernel_size && co[0] <= half_kernel_size ){

        // Compute kernel index 
        // - keeping in mind the central elements are missing (forced to 0)
        //
        
        int kernel_idx = co_to_idx( co+intd2(half_kernel_size, half_kernel_size), intd2(kernel_size, kernel_size) );
        
        if( co[1] == 0 && co[0] == 0 ) {
          kspace[idx] = float_complext(0.0f);
        }
        else {
          if( (co[1]==0 && co[0]>0) || co[1]>0 ) kernel_idx--;
          kspace[idx] = kernels[batch*num_kernel_elements*num_coils + coil*num_kernel_elements + kernel_idx];
        }
      }
      else{
        kspace[idx] = float_complext(0.0f);
      }          
    }
  }

  boost::shared_ptr< cuNDArray<float_complext> > 
  estimate_spirit_kernels( cuNDArray<float_complext> *_kspace, unsigned int kernel_size )
  {
    // Calibration is performed in k-space. 
    // The result is Fourier transformed and returned as image space kernels.
    // The convolution is expressed explicitly as a matrix equation an solved using BLAS/LAPACK.
    //

    if( _kspace == 0x0 ){
      throw std::runtime_error("estimate_spirit_kernels: 0x0 input array");
    }
    
    if( _kspace->get_number_of_dimensions() != 3 ) {
      throw std::runtime_error("estimate_spirit_kernels: Only 2D spirit is supported currently");
    }

    if( (kernel_size%2) == 0 ) {
      throw std::runtime_error("estimate_spirit_kernels: The kernel size should be odd");
    }


    // Normalize input array to an average intensity of one per element
    //
    std::vector<size_t> old_dims = *_kspace->get_dimensions();
    std::vector<size_t> dims= old_dims;
    /*dims[0] /= 2;
    dims[1] /= 2;*/
    //dims[0]=36;
    //dims[1]=36;
    //cuNDArray<float_complext> kspace(_kspace);

    vector_td<size_t,2> offset((old_dims[0]-dims[0])/2,(old_dims[1]-dims[1])/2);
    cuNDArray<float_complext> kspace = crop<float_complext,2>(offset,from_std_vector<size_t,2>(dims),*_kspace);
    float sum = nrm2(&kspace);
    float_complext in_max = kspace[amax(&kspace)];
    kspace /= (float(kspace.get_number_of_elements())/sum);
    unsigned int num_coils = kspace.get_size(kspace.get_number_of_dimensions()-1);
    unsigned int elements_per_coil = kspace.get_number_of_elements()/num_coils;
    
    std::vector<size_t> out_dims;
    out_dims.push_back(_kspace->get_size(0)); out_dims.push_back(_kspace->get_size(1));
    out_dims.push_back(num_coils*num_coils);
    
    boost::shared_ptr< cuNDArray<float_complext> > kernel_images
      ( new cuNDArray<float_complext>(&out_dims) );

    // Clear to ones in case we terminate early
    //

    fill(kernel_images.get(), float_complext(1.0f/num_coils));

    // Form m x n system matrix A
    //

    unsigned int m = elements_per_coil;
    unsigned int n = num_coils*(kernel_size*kernel_size-1);

    std::vector<size_t> A_dims; A_dims.push_back(m); A_dims.push_back(n);    
    cuNDArray<float_complext> A(&A_dims); clear(&A);

    // Fill system matrix
    //

    dim3 blockDim; dim3 gridDim;
    setup_grid( kspace.get_number_of_elements(), &blockDim, &gridDim );
    
    compute_system_matrix_kernel<<< gridDim, blockDim >>>
      ( intd2(kspace.get_size(0), kspace.get_size(1)), num_coils, kernel_size,
        kspace.get_data_ptr(), A.get_data_ptr() );

    CHECK_FOR_CUDA_ERROR();    

    /*
    static int counter = 0;
    char filename[256];
    sprintf((char*)filename, "_A_%d.cplx", counter);
    write_nd_array<float_complext>( A.to_host().get(), filename );
    counter++;
    */

    // Compute A^H A
    //

    cublasStatus_t stat;
    cublasHandle_t handle = *CUBLASContextProvider::instance()->getCublasHandle();

    std::vector<size_t> AHA_dims(2,n);
    cuNDArray<float_complext> AHA(&AHA_dims);

    // Initialize AHA to identity (Tikhonov regularization)
    //

    float_complext one(1.0f);
    clear(&AHA);
    for( unsigned int i=0; i<n; i++ ){
      cudaMemcpy( AHA.get_data_ptr()+i*n+i, &one, sizeof(float_complext), cudaMemcpyHostToDevice );
    }
    CHECK_FOR_CUDA_ERROR();

    float_complext alpha(1.0f);
    //float_complext beta(0.1f*in_max); // Tikhonov regularization weight
    float_complext beta(0.0f); // Tikhonov regularization weight
    
    stat = cublasCgemm( handle, CUBLAS_OP_C, CUBLAS_OP_N,
                        n,n,m,
                        (cuFloatComplex*) &alpha,
                        (cuFloatComplex*) A.get_data_ptr(), m,
                        (cuFloatComplex*) A.get_data_ptr(), m,
                        (cuFloatComplex*) &beta, 
                        (cuFloatComplex*) AHA.get_data_ptr(), n );
    
    if (stat != CUBLAS_STATUS_SUCCESS) {
      std::cerr << "CUBLAS error code " << stat << std::endl;
      throw std::runtime_error("estimate_spirit_kernels: CUBLAS error computing A^HA");
    }

    /*
    static int counter = 0;
    char filename[256];
    sprintf((char*)filename, "_AHA_%d.cplx", counter);
    write_nd_array<float_complext>( AHA.to_host().get(), filename );
    counter++;
    */

    // Multiply A^H with each coil image (to form the rhs)
    //

    std::vector<size_t> rhs_dims; rhs_dims.push_back(n); rhs_dims.push_back(num_coils);    
    cuNDArray<float_complext> rhs(&rhs_dims); clear(&rhs);

    beta = float_complext(0.0f);

    stat = cublasCgemm( handle, CUBLAS_OP_C, CUBLAS_OP_N,
                        n, num_coils, m,
                        (cuFloatComplex*) &alpha,
                        (cuFloatComplex*) A.get_data_ptr(), m,
                        (cuFloatComplex*) kspace.get_data_ptr(), m,
                        (cuFloatComplex*) &beta, 
                        (cuFloatComplex*) rhs.get_data_ptr(), n );
    
    if (stat != CUBLAS_STATUS_SUCCESS) {
      std::cerr << "CUBLAS error code " << stat << std::endl;
      throw std::runtime_error("estimate_spirit_kernels: CUBLAS error computing rhs");
    }
    
    /*
    static int counter = 0;
    char filename[256];
    sprintf((char*)filename, "_rhs_%d.cplx", counter);
    write_nd_array<float_complext>( rhs.to_host().get(), filename );
    counter++;
    */



    //CGELS is used rather than a more conventional solver as it is part of CULA free.
    /*
    culaStatus s = culaDeviceCgels( 'N', n, n, num_coils,
                                 (culaDeviceFloatComplex*)AHA.get_data_ptr(), n,
                                 (culaDeviceFloatComplex*)rhs.get_data_ptr(), n);
    */
    {
      //It actually turns out to be faster to do this inversion on the CPU. Problem is probably too small for GPU to make sense
      //GPUTimer cpu_invert_time("CPU Inversion time");
      ht_grappa_solve_spd_system(AHA, rhs);
    }


    /*
    if( s != culaNoError ) {
      if( s == 8 ){
        std::cerr << "CULA error code " << s << ": " << culaGetStatusString(s) << std::endl;
        std::cerr << "Assuming that the buffer is not yet filled and return ones" << std::endl;
        return kernel_images;
      }
      std::cerr << "CULA error code " << s << ": " << culaGetStatusString(s) << std::endl;
      culaInfo i = culaGetErrorInfo();
      char buf[2048];
      culaGetErrorInfoString(s, i, buf, sizeof(buf));
      printf("Error %d: %s\n", (int)i, buf);      
      throw std::runtime_error("estimate_spirit_kernels: CULA error computing 'getrs'");
    }
    */

    //CULA will sometime return NaN without an explicit error. This code tests for NaNs and returns if found.
    float nan_test = nrm2(&rhs);
    if (nan_test != nan_test) return kernel_images;

    // Fill k-spaces with the computed kernels at the center
    //

    setup_grid( kernel_images->get_number_of_elements(), &blockDim, &gridDim );
    
    write_convolution_masks_kernel<<< gridDim, blockDim >>>
      ( intd2(kernel_images->get_size(0), kernel_images->get_size(1)), num_coils, kernel_size,
        rhs.get_data_ptr(), kernel_images->get_data_ptr() );
    
    CHECK_FOR_CUDA_ERROR();

    // Batch FFT into image space
    //
    A.clear();
    AHA.clear();
    rhs.clear();

    std::vector<size_t> dims_to_xform;
    dims_to_xform.push_back(0); dims_to_xform.push_back(1);    
    cuNDFFT<float>::instance()->ifft( kernel_images.get(), &dims_to_xform, false );
    
    /*
    static int counter = 0;
    char filename[256];
    sprintf((char*)filename, "_kernels_%d.cplx", counter);
    write_nd_array<float_complext>( kernel_images->to_host().get(), filename );
    counter++;
    */

    return kernel_images;
  }
}
