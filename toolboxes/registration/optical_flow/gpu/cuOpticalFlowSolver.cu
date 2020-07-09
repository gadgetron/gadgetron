#include "cuOpticalFlowSolver.h"
#include "vector_td_utilities.h"
#include "check_CUDA.h"

#include <stdexcept>

namespace Gadgetron{

  //
  // Kernel prototype declarations
  //

  template<class REAL, unsigned int D> __global__ 
  void spatial_grad_kernel(const REAL*, const REAL*,REAL*,typename uint64d<D>::Type,unsigned int,unsigned int);

  template<class REAL, unsigned int D> __global__ 
  void temporal_grad_kernel(const REAL*, const REAL*,REAL*,typename uint64d<D>::Type,unsigned int,unsigned int);

  // There is some issue about Cuda defining min/max incompatibly...
  //

  template <class T> __host__ __device__ const T& _cuOF_max (const T& a, const T& b) {
    return (a<b)?b:a;
  }

  template <class T> __host__ __device__ const T& _cuOF_min (const T& a, const T& b) {
    return (a>b)?b:a;
  }

  template<class T, unsigned int D> void
  cuOpticalFlowSolver<T,D>::setup_grid( dim3 *blockDim, dim3* gridDim, 
					   unsigned int number_of_elements, 
					   unsigned int num_batches, 
					   bool use_2d_blocks, 
					   unsigned int num_unknowns )
  {
    int device;
    cudaDeviceProp deviceProp; 
  
    if( cudaGetDevice( &device ) != cudaSuccess) {
      throw std::runtime_error("cuOpticalFlowSolver::setup_grid(): unable to determine current device");
    }
    
    if( cudaGetDeviceProperties( &deviceProp, device ) != cudaSuccess) {
      throw std::runtime_error("cuOpticalFlowSolver::setup_grid(): unable to query current device");
    }
    
    int max_blockdim = deviceProp.maxThreadsDim[0];
    int max_griddim  = deviceProp.maxGridSize[0];
    int warp_size    = deviceProp.warpSize;
    
    // For small arrays we keep the block dimension fairly small
    if( use_2d_blocks )
      *blockDim = dim3(((256/num_unknowns)/warp_size)*warp_size, num_unknowns);
    else
      *blockDim = dim3(256);
  
    *gridDim = dim3((number_of_elements+(blockDim->x*blockDim->y)-1)/(blockDim->x*blockDim->y), num_batches);

    // Extend block/grid dimensions for large arrays
    if( gridDim->x > max_griddim ){
      if( use_2d_blocks )
        blockDim->x = ((max_blockdim/num_unknowns)/warp_size)*warp_size;
      else
        blockDim->x = max_blockdim;
    
      gridDim->x = (number_of_elements+(blockDim->x*blockDim->y)-1)/(blockDim->x*blockDim->y);
    }

    if( gridDim->x > max_griddim ){
      gridDim->x = ((unsigned int)std::sqrt((T)number_of_elements)+(blockDim->x*blockDim->y)-1)/(blockDim->x*blockDim->y);
      gridDim->y *= ((number_of_elements+(blockDim->x*blockDim->y)*gridDim->x-1)/((blockDim->x*blockDim->y)*gridDim->x));
    }
   
    if( gridDim->x > max_griddim || gridDim->y > max_griddim ){      
      throw std::runtime_error("cuOpticalFlowSolver::setup_grid(): maximum grid dimensions exceeded");
    }
  }
  
  template<class T, unsigned int D> void
  cuOpticalFlowSolver<T,D>::core_grad_spatial( T *fixed_image, T *moving_image, T *gradient_image, 
						  typename uint64d<D>::Type matrix_size_moving, 
						  size_t number_of_batches_fixed, 
						  size_t number_of_batches_moving )
  {        
    unsigned int number_of_elements = prod(matrix_size_moving);
    dim3 blockDim; dim3 gridDim;

    setup_grid( &blockDim, &gridDim, number_of_elements, _cuOF_max(number_of_batches_moving, number_of_batches_fixed)*D );
    
    // Invoke kernel (spatial partial derivatives)
    spatial_grad_kernel<T,D><<< gridDim, blockDim >>>
      ( fixed_image, moving_image, gradient_image, matrix_size_moving, number_of_batches_fixed, number_of_batches_moving );
    
    CHECK_FOR_CUDA_ERROR();
  }
  
  template<class T, unsigned int D> void
  cuOpticalFlowSolver<T,D>::core_grad_temporal( T *fixed_image, T *moving_image, T *gradient_image, 
						   typename uint64d<D>::Type matrix_size_moving, 
						   size_t number_of_batches_fixed, 
						   size_t number_of_batches_moving )
  {        
    unsigned int number_of_elements = prod(matrix_size_moving);
    dim3 blockDim; dim3 gridDim;
    
    setup_grid( &blockDim, &gridDim, number_of_elements, _cuOF_max(number_of_batches_moving, number_of_batches_fixed) );
    
    // Invoke kernel (temporal partial derivative)
    temporal_grad_kernel<T,D><<< gridDim, blockDim >>>
      ( fixed_image, moving_image, gradient_image,
        matrix_size_moving, number_of_batches_fixed, number_of_batches_moving );
    
    CHECK_FOR_CUDA_ERROR();
  }
  
  // Helpers
  //

  template<unsigned int D> __device__ 
  typename uint64d<D>::Type compute_stride( unsigned int dim )
  {
    typename uint64d<D>::Type res;
  
    for( unsigned int d=0; d<D; d++ ){
      res.vec[d] = (d==dim) ? 1 : 0;
    }
    return res;
  }

  template<unsigned int D> __device__ 
  bool is_border_pixel_in_stride_dim_before( unsigned int dim, typename uint64d<D>::Type co, typename uint64d<D>::Type dims )
  {
    if( co.vec[dim] == 0 )
      return true;
    else
      return false;
  }

  template<unsigned int D> __device__ 
  bool is_border_pixel_in_stride_dim_after( unsigned int dim, typename uint64d<D>::Type co, typename uint64d<D>::Type dims )
  {
    if( co.vec[dim] == (dims.vec[dim]-1) )
      return true;
    else
      return false;
  }

  // Spatial partial derivatives
  //

  template<class REAL, unsigned int D> __global__ void
  spatial_grad_kernel( const REAL * __restrict__ fixed_image, const REAL * __restrict__ moving_image, REAL * __restrict__ gradient_image,
                       typename uint64d<D>::Type matrix_size, 
                       unsigned int num_batches_fixed, unsigned int num_batches_moving )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

    // Number of elements per partial derivate
    const unsigned int num_elements_per_batch = prod(matrix_size);
    const unsigned int num_elements_per_pdev_fixed = num_elements_per_batch*num_batches_fixed;
    const unsigned int num_elements_per_pdev_moving = num_elements_per_batch*num_batches_moving;

    // Total number of elements for all partial derivatives
    const unsigned int num_elements_total = _cuOF_max(num_elements_per_pdev_fixed, num_elements_per_pdev_moving)*D;
  
    if( idx < num_elements_total ){
    
      // The (minimum) index in the slowest varying output dimension determines which partial derivative to compute 
      const unsigned int stride_dim_fixed = idx/(num_elements_per_pdev_fixed);
      const unsigned int stride_dim_moving = idx/(num_elements_per_pdev_moving);
      const unsigned int stride_dim = _cuOF_min(stride_dim_fixed, stride_dim_moving);

      // Local index to the partial derivative
      const unsigned int idx_in_pdev_fixed = idx-stride_dim_fixed*num_elements_per_pdev_fixed;
      const unsigned int idx_in_pdev_moving = idx-stride_dim_moving*num_elements_per_pdev_moving;

      // Batch idx (second slowest varying dimension)   
      const unsigned int batch_idx_fixed = idx_in_pdev_fixed/num_elements_per_batch;
      const unsigned int batch_idx_moving = idx_in_pdev_moving/num_elements_per_batch;

      // Local index to the batch (should be identical for the fixed/moving image)
      const size_t idx_in_batch = idx_in_pdev_moving-batch_idx_moving*num_elements_per_batch;

      // Local co to the image
      const typename uint64d<D>::Type co = idx_to_co( idx_in_batch, matrix_size );
 
      REAL res;
      unsigned int count = 0;

      //
      // Find partial derivatives using central differences
      //
    
      typename uint64d<D>::Type stride = compute_stride<D>(stride_dim);
    
      const unsigned int base_idx_moving = batch_idx_moving*num_elements_per_batch;
      const unsigned int base_idx_fixed = batch_idx_fixed*num_elements_per_batch;

      unsigned int stride_base_idx, fixed_idx, moving_idx;
     
      // Neighbor "plus stride" side
      if( !is_border_pixel_in_stride_dim_after<D>( stride_dim, co, matrix_size )){
        stride_base_idx = co_to_idx(co+stride, matrix_size);
        count++;
      }
      else{
        stride_base_idx = idx_in_batch;
      }
    
      fixed_idx = stride_base_idx+base_idx_fixed;
      moving_idx = stride_base_idx+base_idx_moving;
    
      res = (fixed_image[fixed_idx]+moving_image[moving_idx])*REAL(0.5);

      // Neighbor "minus stride" side
      if( !is_border_pixel_in_stride_dim_before<D>( stride_dim, co, matrix_size )){
        stride_base_idx = co_to_idx(co-stride, matrix_size);
        count++;
      }
      else{
        stride_base_idx = co_to_idx(co, matrix_size);
      }
    
      fixed_idx = stride_base_idx+base_idx_fixed;
      moving_idx = stride_base_idx+base_idx_moving;
    
      res -= (fixed_image[fixed_idx]+moving_image[moving_idx])*REAL(0.5);

      if( count == 2 ) // Both neighbors exist
        res /= REAL(2);

      // Output result
      //
    
      gradient_image[idx] = res;
    }
  }

  // Temporal partial derivatives
  //

  template<class REAL, unsigned int D> __global__ void
  temporal_grad_kernel( const REAL * __restrict__ fixed_image, const REAL * __restrict__ moving_image, REAL * __restrict__ gradient_image,
                        typename uint64d<D>::Type matrix_size, 
                        unsigned int num_batches_fixed, unsigned int num_batches_moving )
  { 
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

    // Number of elements per partial derivate
    const unsigned int num_elements_per_batch = prod(matrix_size);
    const unsigned int num_elements_per_pdev_fixed = num_elements_per_batch*num_batches_fixed;
    const unsigned int num_elements_per_pdev_moving = num_elements_per_batch*num_batches_moving;

    // Total number of elements for all partial derivatives
    const unsigned int num_elements_total = _cuOF_max(num_elements_per_pdev_fixed, num_elements_per_pdev_moving);
  
    if( idx < num_elements_total ){
    
      const unsigned int stride_dim_fixed = idx/(num_elements_per_pdev_fixed);
      const unsigned int stride_dim_moving = idx/(num_elements_per_pdev_moving);

      // Local index to the partial derivative
      const unsigned int idx_in_pdev_fixed = idx-stride_dim_fixed*num_elements_per_pdev_fixed;
      const unsigned int idx_in_pdev_moving = idx-stride_dim_moving*num_elements_per_pdev_moving;

      // Batch idx (second slowest varying dimension)   
      const unsigned int batch_idx_fixed = idx_in_pdev_fixed/num_elements_per_batch;
      const unsigned int batch_idx_moving = idx_in_pdev_moving/num_elements_per_batch;

      // Local index to the batch (should be identical for the fixed/moving image)
      const unsigned int idx_in_batch = idx_in_pdev_moving-batch_idx_moving*num_elements_per_batch;

      const unsigned int base_idx_fixed = batch_idx_fixed*num_elements_per_batch;
      const unsigned int base_idx_moving = batch_idx_moving*num_elements_per_batch;
    
      // Ctr pixel
      const unsigned int fixed_idx = idx_in_batch+base_idx_fixed;
      const unsigned int moving_idx = idx_in_batch+base_idx_moving;
    
      const REAL res = moving_image[moving_idx]-fixed_image[fixed_idx];
    
      // Output result
      //
    
      gradient_image[idx] = res;        
    }    
  }

  // 
  // Template instantiation
  //

  template class EXPORTGPUREG cuOpticalFlowSolver<float,1>;
  template class EXPORTGPUREG cuOpticalFlowSolver<float,2>;
  template class EXPORTGPUREG cuOpticalFlowSolver<float,3>;
  template class EXPORTGPUREG cuOpticalFlowSolver<float,4>;

  template class EXPORTGPUREG cuOpticalFlowSolver<double,1>;
  template class EXPORTGPUREG cuOpticalFlowSolver<double,2>;
  template class EXPORTGPUREG cuOpticalFlowSolver<double,3>;
  template class EXPORTGPUREG cuOpticalFlowSolver<double,4>;
}
