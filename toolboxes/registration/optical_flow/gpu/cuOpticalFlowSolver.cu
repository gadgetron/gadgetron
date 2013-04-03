#include "cuOpticalFlowSolver.h"
#include "ndarray_vector_td_utilities.h"
#include "check_CUDA.h"

//
// Kernel prototype declarations
//

template<class REAL, unsigned int D> __global__ 
void spatial_grad_kernel(REAL*,REAL*,REAL*,typename uintd<D>::Type,unsigned int,unsigned int);

template<class REAL, unsigned int D> __global__ 
void temporal_grad_kernel(REAL*,REAL*,REAL*,typename uintd<D>::Type,unsigned int,unsigned int);

//
// Implementation
//

template<class REAL, unsigned int D> bool
cuOpticalFlowSolver<REAL,D>::normalize( cuNDArray<REAL> *image )
{
  REAL scale = cuNDA_normalize<REAL>(image, REAL(1));
  if( scale == REAL(0) )
    return false;
  else
    return true;
}

template<class REAL, unsigned int D> boost::shared_ptr< cuNDArray<REAL> > 
cuOpticalFlowSolver<REAL,D>::downsample( cuNDArray<REAL> *image )
{
  return cuNDA_downsample<REAL,D>(image);
}

template<class REAL, unsigned int D> boost::shared_ptr< cuNDArray<REAL> > 
cuOpticalFlowSolver<REAL,D>::upsample( cuNDArray<REAL> *displacements )
{
  boost::shared_ptr< cuNDArray<REAL> > disp_hires = cuNDA_upsample_lin<REAL,D>(displacements);
  
  // After upsampling the vector field should be multiplied by 2 to match the higher resolution
  //

  if( !cuNDA_scal<REAL>( REAL(2.0), disp_hires.get() )){
    this->solver_error( "cuOpticalFlowSolver::upsample : failed to rescale array" );
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  return disp_hires;
} 

template<class REAL, unsigned int D>  bool 
cuOpticalFlowSolver<REAL,D>::compute( cuNDArray<REAL> *fixed_image, 
				      cuNDArray<REAL> *moving_image, 
				      cuNDArray<REAL> *stencil_image, 
				      boost::shared_ptr< cuNDArray<REAL> > &result_in_out )
{
  // Test the validity of the input images
  //
  
  if( !fixed_image || !moving_image ){
    this->solver_error( "cuOpticalFlowSolver::compute : NULL input images provided" );
    return false;
  }

  if( prod(vector_to_uintd<D>(*fixed_image->get_dimensions().get())) != 
      prod(vector_to_uintd<D>(*moving_image->get_dimensions().get())) ){
    this->solver_error( "cuOpticalFlowSolver::compute : core image dimensions (excluding batches) mismatch" );
    return false;
  }

  if( stencil_image && 
      prod(vector_to_uintd<D>(*fixed_image->get_dimensions().get())) != 
      prod(vector_to_uintd<D>(*stencil_image->get_dimensions().get())) ){
    this->solver_error( "cuOpticalFlowSolver::compute : stencil image dimensions mismatch fixed/moving image dimensions" );
    return false;
  }
  
  if( result_in_out.get() && 
      !( result_in_out->get_number_of_dimensions() > D ||
	 result_in_out->get_size(result_in_out->get_number_of_dimensions()-1) == D )){
    this->solver_error( "cuOpticalFlowSolver::compute : input displacements dimensionality mismatch" );
    return false;
  }
  
  // If an approximate displacement field is provided it is used to resample the moving image
  //
  
  boost::shared_ptr< cuNDArray<REAL> > _def_moving_image;
  cuNDArray<REAL> *def_moving_image = 0x0;

  if( result_in_out.get() ){ 

    // Apply the input deformation
    //
    
    _def_moving_image = deform( moving_image, result_in_out );
    def_moving_image = _def_moving_image.get();
  }
  else{
    
    // The is no input deformation to apply
    //

    def_moving_image = moving_image;
  }
  
  // Compute gradient image
  //

  boost::shared_ptr< cuNDArray<REAL> > grad_image = grad( fixed_image, def_moving_image );
  if( !grad_image.get() ){
    this->solver_error( "cuOpticalFlowSolver::compute : gradient image computation failed" );
    return false;
  }

  // The deformed image is no longer needed
  //

  _def_moving_image.reset(); def_moving_image = 0x0;
  
  // Invoke core solver (e.g. Horn-Schunk, Cornelius-Kanade, ...)
  //

  boost::shared_ptr< cuNDArray<REAL> > displacements = core_solver( grad_image.get(), stencil_image );
    
  if( !displacements.get() ){
    this->solver_error( "cuOpticalFlowSolver::compute : core optical flow solver failed" );
    return false;
  }
  
  // If an input vector field was provided then our result should be added element-wise
  // 
  
  if( result_in_out.get() ){
    cuNDA_axpy( REAL(1), displacements.get(), result_in_out.get() );
  }
  else{    
    result_in_out = displacements;
  }
  
  return true;  
}
    
template<class REAL, unsigned int D> boost::shared_ptr< cuNDArray<REAL> > 
cuOpticalFlowSolver<REAL,D>::grad( cuNDArray<REAL> *fixed_image, cuNDArray<REAL> *moving_image )
{
  // Sanity checks
  //
  
  if( !fixed_image || !moving_image ){
    this->solver_error( "cuOpticalFlowSolver::grad : NULL input" );
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  if( !((moving_image->get_number_of_elements() % fixed_image->get_number_of_elements()) == 0 ||
	(fixed_image->get_number_of_elements() % moving_image->get_number_of_elements()) == 0 )){
    this->solver_error( "cuOpticalFlowSolver::grad : fixed/moving image dimensions mismatch" );
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Determine dimension size of the gradient field:
  // D spatial dimensions plus one temporal dimension
  //
  
  std::vector<unsigned int> grad_dims;

  (fixed_image->get_number_of_elements()<moving_image->get_number_of_elements() )
    ? grad_dims = *(moving_image->get_dimensions().get())
    : grad_dims = *(fixed_image->get_dimensions().get());
  
  grad_dims.push_back(D+1); 
  
  boost::shared_ptr< cuNDArray<REAL> > grad_image( new cuNDArray<REAL>());
  if( !grad_image->create( &grad_dims )) {
    this->solver_error( "cuOpticalFlowSolver::grad : Unable to allocate storage for gradient image" );
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Setup for the spatial partial derivatives
  //
  
  typename uintd<D>::Type matrix_size_fixed = vector_to_uintd<D>( *fixed_image->get_dimensions() );
  typename uintd<D>::Type matrix_size_moving = vector_to_uintd<D>( *moving_image->get_dimensions() );

  if( matrix_size_fixed != matrix_size_moving ){
    this->solver_error( "cuOpticalFlowSolver::grad : fixed and moving image dimensions (ignoring batch dims) mismatch" );
    return boost::shared_ptr< cuNDArray<REAL> >();
  }
  
  // Ignoring the batch dimensions the fixed and moving images have the same number of elements
  //
  
  unsigned int number_of_elements = prod(matrix_size_moving);
  
  unsigned int number_of_batches_fixed = 1;
  unsigned int number_of_batches_moving = 1;

  for( unsigned int d=D; d<fixed_image->get_number_of_dimensions(); d++ ){
    number_of_batches_fixed *= fixed_image->get_size(d);
  }
  
  for( unsigned int d=D; d<moving_image->get_number_of_dimensions(); d++ ){
    number_of_batches_moving *= moving_image->get_size(d);
  }
  
  // Determine grid configuration (spatial partial derivatives)
  //

  dim3 blockDim; dim3 gridDim;
  if( !setup_grid( &blockDim, &gridDim, number_of_elements, 
		   max(number_of_batches_moving, number_of_batches_fixed)*D )){
    this->solver_error( "cuOpticalFlowSolver::grad : Unable determine grid dimensions (1)" );
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Invoke kernel (spatial partial derivatives)
  spatial_grad_kernel<REAL,D><<< gridDim, blockDim >>>
    ( fixed_image->get_data_ptr(), moving_image->get_data_ptr(), grad_image->get_data_ptr(), 
      matrix_size_moving, number_of_batches_fixed, number_of_batches_moving );
  
  CHECK_FOR_CUDA_ERROR();
  
  // Determine grid configuration (temporal partial derivative)
  //
  
  if( !setup_grid( &blockDim, &gridDim, number_of_elements, 
		   max(number_of_batches_moving, number_of_batches_fixed)*1 )){
    this->solver_error( "cuOpticalFlowSolver::grad : Unable determine grid dimensions (2)" );
    return boost::shared_ptr< cuNDArray<REAL> >();
  }

  // Invoke kernel (temporal partial derivative)
  temporal_grad_kernel<REAL,D><<< gridDim, blockDim >>>
    ( fixed_image->get_data_ptr(), moving_image->get_data_ptr(), 
      grad_image->get_data_ptr()+number_of_elements*max(number_of_batches_moving, number_of_batches_fixed)*D, 
      matrix_size_moving, number_of_batches_fixed, number_of_batches_moving );
  
  CHECK_FOR_CUDA_ERROR();
  
  // We are done...
  //

  return grad_image;
}

template<class REAL, unsigned int D> bool 
cuOpticalFlowSolver<REAL,D>::setup_grid( dim3 *blockDim, dim3* gridDim, unsigned int number_of_elements, unsigned int num_batches, bool use_2d_blocks, unsigned int num_unknowns )
{
  int device;
  cudaDeviceProp deviceProp; 
  
  if( cudaGetDevice( &device ) != cudaSuccess) {
    std::cout << std::endl << "Error: unable to determine current device." << std::endl ;
    return false;
  }

  if( cudaGetDeviceProperties( &deviceProp, device ) != cudaSuccess) {
    std::cout << std::endl << "Error: unable to query device properties." << std::endl ;
    return false;
  }
  
  int max_blockdim = deviceProp.maxThreadsDim[0];
  int max_griddim  = deviceProp.maxGridSize[0];
  int warp_size    = deviceProp.warpSize;
  
  // For small arrays we keep the block dimension fairly small
  if( use_2d_blocks )
    *blockDim = dim3(((256/num_unknowns)/warp_size)*warp_size,num_unknowns);
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
    gridDim->x = ((unsigned int)sqrt((float)number_of_elements)+(blockDim->x*blockDim->y)-1)/(blockDim->x*blockDim->y);
    gridDim->y *= ((number_of_elements+(blockDim->x*blockDim->y)*gridDim->x-1)/((blockDim->x*blockDim->y)*gridDim->x));
  }
   
  if( gridDim->x > max_griddim || gridDim->y > max_griddim )
    return false;
  else 
    return true;
}

//
// Kernels
//

// Helpers
//

template<unsigned int D> __device__ 
typename uintd<D>::Type compute_stride( unsigned int dim )
{
  typename uintd<D>::Type res;
  
  for( unsigned int d=0; d<D; d++ ){
    res.vec[d] = (d==dim) ? 1 : 0;
  }
  return res;
}

template<unsigned int D> __device__ 
bool is_border_pixel_in_stride_dim_before( unsigned int dim, typename uintd<D>::Type co, typename uintd<D>::Type dims )
{
  if( co.vec[dim] == 0 )
    return true;
  else
    return false;
}

template<unsigned int D> __device__ 
bool is_border_pixel_in_stride_dim_after( unsigned int dim, typename uintd<D>::Type co, typename uintd<D>::Type dims )
{
  if( co.vec[dim] == (dims.vec[dim]-1) )
    return true;
  else
    return false;
}

// Spatial partial derivatives
//

template<class REAL, unsigned int D> __global__ void
spatial_grad_kernel( REAL *fixed_image, REAL *moving_image, REAL *gradient_image, 
		     typename uintd<D>::Type matrix_size, 
		     unsigned int num_batches_fixed, unsigned int num_batches_moving )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  // Number of elements per partial derivate
  const unsigned int num_elements_per_batch = prod(matrix_size);
  const unsigned int num_elements_per_pdev_fixed = num_elements_per_batch*num_batches_fixed;
  const unsigned int num_elements_per_pdev_moving = num_elements_per_batch*num_batches_moving;

  // Total number of elements for all partial derivatives
  const unsigned int num_elements_total = max(num_elements_per_pdev_fixed, num_elements_per_pdev_moving)*D;
  
  if( idx < num_elements_total ){
    
    // The (minimum) index in the slowest varying output dimension determines which partial derivative to compute 
    const unsigned int stride_dim_fixed = idx/(num_elements_per_pdev_fixed);
    const unsigned int stride_dim_moving = idx/(num_elements_per_pdev_moving);
    const unsigned int stride_dim = min(stride_dim_fixed, stride_dim_moving);

    // Local index to the partial derivative
    const unsigned int idx_in_pdev_fixed = idx-stride_dim_fixed*num_elements_per_pdev_fixed;
    const unsigned int idx_in_pdev_moving = idx-stride_dim_moving*num_elements_per_pdev_moving;

    // Batch idx (second slowest varying dimension)   
    const unsigned int batch_idx_fixed = idx_in_pdev_fixed/num_elements_per_batch;
    const unsigned int batch_idx_moving = idx_in_pdev_moving/num_elements_per_batch;

    // Local index to the batch (should be identical for the fixed/moving image)
    const unsigned int idx_in_batch = idx_in_pdev_moving-batch_idx_moving*num_elements_per_batch;

    // Local co to the image
    const typename uintd<D>::Type co = idx_to_co<D>( idx_in_batch, matrix_size );
 
    REAL res;
    unsigned int count = 0;

    //
    // Find partial derivatives       
    // Use central differences
    //
    
    typename uintd<D>::Type stride = compute_stride<D>(stride_dim);
    
    const unsigned int base_idx_moving = batch_idx_moving*num_elements_per_batch;
    const unsigned int base_idx_fixed = batch_idx_fixed*num_elements_per_batch;

    unsigned int stride_base_idx, fixed_idx, moving_idx;
     
    // Neighbor "plus stride" side
    if( !is_border_pixel_in_stride_dim_after<D>( stride_dim, co, matrix_size )){
      stride_base_idx = co_to_idx<D>(co+stride, matrix_size);
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
      stride_base_idx = co_to_idx<D>(co-stride, matrix_size);
      count++;
    }
    else{
      stride_base_idx = co_to_idx<D>(co, matrix_size);
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
temporal_grad_kernel( REAL *fixed_image, REAL *moving_image, REAL *gradient_image, 
		      typename uintd<D>::Type matrix_size, 
		      unsigned int num_batches_fixed, unsigned int num_batches_moving )
{ 
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  // Number of elements per partial derivate
  const unsigned int num_elements_per_batch = prod(matrix_size);
  const unsigned int num_elements_per_pdev_fixed = num_elements_per_batch*num_batches_fixed;
  const unsigned int num_elements_per_pdev_moving = num_elements_per_batch*num_batches_moving;

  // Total number of elements for all partial derivatives
  const unsigned int num_elements_total = max(num_elements_per_pdev_fixed, num_elements_per_pdev_moving)*D;
  
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

