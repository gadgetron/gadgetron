#include "hoOpticalFlowSolver.h"
#include "vector_td_utilities.h"

#include <algorithm>

namespace Gadgetron{

  // Helpers
  //

  template<unsigned int D> inline typename uint64d<D>::Type 
  compute_stride( size_t dim )
  {
    typename uint64d<D>::Type res;
  
    for( size_t d=0; d<D; d++ ){
      res.vec[d] = (d==dim) ? 1 : 0;
    }
    return res;
  }

  template<unsigned int D> inline bool 
  is_border_pixel_in_stride_dim_before( size_t dim, typename uint64d<D>::Type co, typename uint64d<D>::Type dims )
  {
    if( co.vec[dim] == 0 )
      return true;
    else
      return false;
  }

  template<unsigned int D> inline bool 
  is_border_pixel_in_stride_dim_after( size_t dim, typename uint64d<D>::Type co, typename uint64d<D>::Type dims )
  {
    if( co.vec[dim] == (dims.vec[dim]-1) )
      return true;
    else
      return false;
  }
    
  template<class T, unsigned int D> void
  hoOpticalFlowSolver<T,D>::core_grad_spatial( T *fixed_image, T *moving_image, T *gradient_image, 
						  typename uint64d<D>::Type matrix_size, 
						  size_t num_batches_fixed, 
						  size_t num_batches_moving )
  {        
    // Number of elements per partial derivate
    const size_t num_elements_per_batch = prod(matrix_size);
    const size_t num_elements_per_pdev_fixed = num_elements_per_batch*num_batches_fixed;
    const size_t num_elements_per_pdev_moving = num_elements_per_batch*num_batches_moving;

    // Total number of elements for all partial derivatives
    const size_t num_elements_total = std::max(num_elements_per_pdev_fixed, num_elements_per_pdev_moving)*D;
  
    for( size_t idx = 0; idx<num_elements_total; idx++ ){
    
      // The (minimum) index in the slowest varying output dimension determines which partial derivative to compute 
      const size_t stride_dim_fixed = idx/(num_elements_per_pdev_fixed);
      const size_t stride_dim_moving = idx/(num_elements_per_pdev_moving);
      const size_t stride_dim = std::min(stride_dim_fixed, stride_dim_moving);

      // Local index to the partial derivative
      const size_t idx_in_pdev_fixed = idx-stride_dim_fixed*num_elements_per_pdev_fixed;
      const size_t idx_in_pdev_moving = idx-stride_dim_moving*num_elements_per_pdev_moving;

      // Batch idx (second slowest varying dimension)   
      const size_t batch_idx_fixed = idx_in_pdev_fixed/num_elements_per_batch;
      const size_t batch_idx_moving = idx_in_pdev_moving/num_elements_per_batch;

      // Local index to the batch (should be identical for the fixed/moving image)
      const size_t idx_in_batch = idx_in_pdev_moving-batch_idx_moving*num_elements_per_batch;

      // Local co to the image
      const typename uint64d<D>::Type co = idx_to_co( idx_in_batch, matrix_size );
 
      T res;
      size_t count = 0;

      //
      // Find partial derivatives using central differences
      //
    
      const typename uint64d<D>::Type stride = compute_stride<D>(stride_dim);
      const size_t base_idx_moving = batch_idx_moving*num_elements_per_batch;
      const size_t base_idx_fixed = batch_idx_fixed*num_elements_per_batch;

      size_t stride_base_idx, fixed_idx, moving_idx;
     
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
      
      res = (fixed_image[fixed_idx]+moving_image[moving_idx])*T(0.5);

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
    
      res -= (fixed_image[fixed_idx]+moving_image[moving_idx])*T(0.5);

      if( count == 2 ) // Both neighbors exist
        res /= T(2);

      // Output result
      //
    
      gradient_image[idx] = res;
    }
  }
  
  template<class T, unsigned int D> void
  hoOpticalFlowSolver<T,D>::core_grad_temporal( T *fixed_image, T *moving_image, T *gradient_image, 
						   typename uint64d<D>::Type matrix_size, 
						   size_t num_batches_fixed, 
						   size_t num_batches_moving )
  {        
    // Number of elements per partial derivate
    const size_t num_elements_per_batch = prod(matrix_size);
    const size_t num_elements_per_pdev_fixed = num_elements_per_batch*num_batches_fixed;
    const size_t num_elements_per_pdev_moving = num_elements_per_batch*num_batches_moving;

    // Total number of elements for all partial derivatives
    const size_t num_elements_total = std::max(num_elements_per_pdev_fixed, num_elements_per_pdev_moving);
  
    for( size_t idx =0; idx < num_elements_total; idx++ ){
      
      // Local index to the partial derivative
      const size_t stride_dim_fixed = idx/(num_elements_per_pdev_fixed);
      const size_t stride_dim_moving = idx/(num_elements_per_pdev_moving);
      const size_t idx_in_pdev_fixed = idx-stride_dim_fixed*num_elements_per_pdev_fixed;
      const size_t idx_in_pdev_moving = idx-stride_dim_moving*num_elements_per_pdev_moving;

      // Batch idx (second slowest varying dimension)   
      const size_t batch_idx_fixed = idx_in_pdev_fixed/num_elements_per_batch;
      const size_t batch_idx_moving = idx_in_pdev_moving/num_elements_per_batch;

      // Local index to the batch (should be identical for the fixed/moving image)
      const size_t idx_in_batch = idx_in_pdev_moving-batch_idx_moving*num_elements_per_batch;
      const size_t base_idx_fixed = batch_idx_fixed*num_elements_per_batch;
      const size_t base_idx_moving = batch_idx_moving*num_elements_per_batch;
    
      // Ctr pixel
      const size_t fixed_idx = idx_in_batch+base_idx_fixed;
      const size_t moving_idx = idx_in_batch+base_idx_moving;
    
      const T res = moving_image[moving_idx]-fixed_image[fixed_idx];
    
      // Output result
      //
    
      gradient_image[idx] = res;        
    }    
  }
  
  // 
  // Template instantiation
  //

  template class EXPORTCPUREG hoOpticalFlowSolver<float,1>;
  template class EXPORTCPUREG hoOpticalFlowSolver<float,2>;
  template class EXPORTCPUREG hoOpticalFlowSolver<float,3>;
  template class EXPORTCPUREG hoOpticalFlowSolver<float,4>;

  template class EXPORTCPUREG hoOpticalFlowSolver<double,1>;
  template class EXPORTCPUREG hoOpticalFlowSolver<double,2>;
  template class EXPORTCPUREG hoOpticalFlowSolver<double,3>;
  template class EXPORTCPUREG hoOpticalFlowSolver<double,4>;
}
