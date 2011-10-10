#include "b1_map.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "real_utilities.h"
#include "real_utilities_device.h"
#include "check_CUDA.h"
#include "cuNDFFT.h"

#include <iostream>

using namespace std;

const int kernel_width = 7;

template<class REAL, unsigned int D> void
smooth_correlation_matrices( cuNDArray<typename complext<REAL>::Type> *corrm, cuNDArray<typename complext<REAL>::Type> *corrm_smooth );

template<class REAL> __host__ 
boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> > extract_csm( cuNDArray<typename complext<REAL>::Type> *corrm_in, unsigned int number_of_batches, unsigned int number_of_elements );

template<class REAL> __host__ 
void set_phase_reference( cuNDArray<typename complext<REAL>::Type> *csm, unsigned int number_of_batches, unsigned int number_of_elements );

//
// Main method
//

template<class REAL, unsigned int D> boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> >
estimate_b1_map( cuNDArray<typename complext<REAL>::Type> *data_in )
{
  if( data_in->get_number_of_dimensions() < 2 ){
    cout << endl << "estimate_b1_map:: dimensionality mismatch." << endl; 
    return boost::shared_ptr< cuNDArray<typename complext<REAL>::Type > >();
  }

  if( data_in->get_number_of_dimensions()-1 != D ){
    cout << endl << "estimate_b1_map:: dimensionality mismatch." << endl; 
    return boost::shared_ptr< cuNDArray<typename complext<REAL>::Type > >();
  }

  vector<unsigned int> image_dims, dims_to_xform;
  unsigned int pixels_per_coil = 1;
  
  for( unsigned int i=0; i<D; i++ ){
    image_dims.push_back(data_in->get_size(i));
    dims_to_xform.push_back(i);
    pixels_per_coil *= data_in->get_size(i);
  }
  
  unsigned int ncoils = data_in->get_size(D);

  // Make a copy of input data
  cuNDArray<typename complext<REAL>::Type > *_data_out = new cuNDArray<typename complext<REAL>::Type>(*data_in);
  boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> > data_out(_data_out);
  
  // Normalize by the RSS of the coils
  if( !cuNDA_rss_normalize<REAL>( data_out.get(), D ) ){
    cout << endl << "estimate_b1_map:: error in rss_normalize" << endl;
    return boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> >();
  }
  
  // Now calculate the correlation matrices
  boost::shared_ptr<cuNDArray<typename complext<REAL>::Type> > corrm = cuNDA_correlation<typename complext<REAL>::Type>( data_out.get() );
  data_out.reset();
  
  // Smooth (onto copy of corrm)
  cuNDArray<typename complext<REAL>::Type > *_corrm_smooth = new cuNDArray<typename complext<REAL>::Type>();
  _corrm_smooth->create(corrm->get_dimensions().get());
  boost::shared_ptr<cuNDArray<typename complext<REAL>::Type> > corrm_smooth(_corrm_smooth);

  smooth_correlation_matrices<REAL,D>( corrm.get(), corrm_smooth.get() );
  corrm.reset();

  // Get the dominant eigenvector for each correlation matrix.
  boost::shared_ptr<cuNDArray<typename complext<REAL>::Type> > csm = extract_csm<REAL>( corrm_smooth.get(), ncoils, pixels_per_coil );
  corrm_smooth.reset();
  
  // Set phase according to reference (coil 0)
  set_phase_reference<REAL>( csm.get(), ncoils, pixels_per_coil );
  
  return csm;
}

// Smooth correlation matrices by box filter (1D)
template<class REAL> __global__ void
smooth_correlation_matrices_kernel( typename complext<REAL>::Type *corrm, typename complext<REAL>::Type *corrm_smooth, intd<1>::Type image_dims )
{
  const int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const int batch = blockIdx.y;

  const int num_image_elements = prod(image_dims);

  if( idx < num_image_elements ){
    
    const int co = idx;    
    const int x = co;
    
    const int size_x = image_dims.vec[0];
    
    const REAL scale = get_one<REAL>()/((REAL)kernel_width);
    
    typename complext<REAL>::Type result = get_zero<typename complext<REAL>::Type>();
    
    for (int kx = 0; kx < kernel_width; kx++) {
      
      if ((x-(kernel_width>>1)+kx) >= 0 &&
	  (x-(kernel_width>>1)+kx) < size_x)
	{	    
	  int source_offset = 
	    batch*num_image_elements +
	    (x-(kernel_width>>1)+kx);
	  
	  result += corrm[source_offset];
	}
    }
    corrm_smooth[batch*num_image_elements+idx] = scale*result;
  }
}

// Smooth correlation matrices by box filter (2D)
template<class REAL> __global__ void
smooth_correlation_matrices_kernel( typename complext<REAL>::Type *corrm, typename complext<REAL>::Type *corrm_smooth, intd<2>::Type image_dims )
{
  const int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const int batch = blockIdx.y;

  const int num_image_elements = prod(image_dims);

  if( idx < num_image_elements ){
    
    const intd2::Type co = idx_to_co<2>(idx, image_dims);
    
    const int x = co.vec[0];
    const int y = co.vec[1];
    
    const int size_x = image_dims.vec[0];
    const int size_y = image_dims.vec[1];
    
    const int half_width = kernel_width>>1;

    const int yminus = y-half_width;
    const int xminus = x-half_width;
    const int yplus = y+half_width;
    const int xplus = x+half_width;

    const REAL scale = get_one<REAL>()/((REAL)(kernel_width*kernel_width));
    
    typename complext<REAL>::Type result = get_zero<typename complext<REAL>::Type>();
   
    if( (yminus >=0) ){
      if( yplus < size_y ){
	if( xminus >= 0 ){
	  if( xplus < size_x ){

#pragma unroll
	    for (int ky = 0; ky < kernel_width; ky++){
#pragma unroll
	      for (int kx = 0; kx < kernel_width; kx++) {
		
		int cy = yminus+ky;
		int cx = xminus+kx;
		
		int source_offset = batch*num_image_elements + cy*size_x + cx;
		result += corrm[source_offset];
	      }
	    }
	  }
	}
      }
    }
    corrm_smooth[batch*num_image_elements+idx] = scale*result;
  }
}

// Smooth correlation matrices by box filter (3D)
template<class REAL> __global__ void
smooth_correlation_matrices_kernel( typename complext<REAL>::Type *corrm, typename complext<REAL>::Type *corrm_smooth, intd<3>::Type image_dims )
{
  const int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const int batch = blockIdx.y;

  const int num_image_elements = prod(image_dims);

  if( idx < num_image_elements ){
    
    const intd3::Type co = idx_to_co<3>(idx, image_dims);
    
    const int x = co.vec[0];
    const int y = co.vec[1];
    const int z = co.vec[2];
    
    const int size_x = image_dims.vec[0];
    const int size_y = image_dims.vec[1];
    const int size_z = image_dims.vec[2];
    
    const REAL scale = get_one<REAL>()/((REAL)(kernel_width*kernel_width*kernel_width));
    
    typename complext<REAL>::Type result = get_zero<typename complext<REAL>::Type>();
    
    for (int kz = 0; kz < kernel_width; kz++) {
      for (int ky = 0; ky < kernel_width; ky++) {
	for (int kx = 0; kx < kernel_width; kx++) {
	
	  if ((z-(kernel_width>>1)+kz) >= 0 &&
	      (z-(kernel_width>>1)+kz) < size_z &&
	      (y-(kernel_width>>1)+ky) >= 0 &&
	      (y-(kernel_width>>1)+ky) < size_y &&
	      (x-(kernel_width>>1)+kx) >= 0 &&
	      (x-(kernel_width>>1)+kx) < size_x) 
	    {	    
	      int source_offset = 
		batch*num_image_elements +
		(z-(kernel_width>>1)+kz)*size_x*size_y +
		(y-(kernel_width>>1)+ky)*size_x +
		(x-(kernel_width>>1)+kx);
	    
	      result += corrm[source_offset];
	    }
	}
      }
    }
    corrm_smooth[batch*num_image_elements+idx] = scale*result;
  }
}

// Smooth correlation matrices by box filter (3D)
template<class REAL> __global__ void
smooth_correlation_matrices_kernel( typename complext<REAL>::Type *corrm, typename complext<REAL>::Type *corrm_smooth, intd<4>::Type image_dims )
{
  const int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const int batch = blockIdx.y;

  const int num_image_elements = prod(image_dims);

  if( idx < num_image_elements ){
    
    const intd4::Type co = idx_to_co<4>(idx, image_dims);
    
    const int x = co.vec[0];
    const int y = co.vec[1];
    const int z = co.vec[2];
    const int w = co.vec[3];
    
    const int size_x = image_dims.vec[0];
    const int size_y = image_dims.vec[1];
    const int size_z = image_dims.vec[2];    
    const int size_w = image_dims.vec[3];
    
    const REAL scale = get_one<REAL>()/((REAL)(kernel_width*kernel_width*kernel_width*kernel_width));
    
    typename complext<REAL>::Type result = get_zero<typename complext<REAL>::Type>();
    
    for (int kw = 0; kw < kernel_width; kw++) {
      for (int kz = 0; kz < kernel_width; kz++) {
	for (int ky = 0; ky < kernel_width; ky++) {
	  for (int kx = 0; kx < kernel_width; kx++) {
	
	    if ((w-(kernel_width>>1)+kw) >= 0 &&
		(w-(kernel_width>>1)+kw) < size_w &&
		(z-(kernel_width>>1)+kz) >= 0 &&
		(z-(kernel_width>>1)+kz) < size_z &&
		(y-(kernel_width>>1)+ky) >= 0 &&
		(y-(kernel_width>>1)+ky) < size_y &&
		(x-(kernel_width>>1)+kx) >= 0 &&
		(x-(kernel_width>>1)+kx) < size_x) 
	      {	    
		int source_offset = 
		  batch*num_image_elements +
		  (w-(kernel_width>>1)+kw)*size_x*size_y*size_z +
		  (z-(kernel_width>>1)+kz)*size_x*size_y +
		  (y-(kernel_width>>1)+ky)*size_x +
		  (x-(kernel_width>>1)+kx);
	    
		result += corrm[source_offset];
	      }
	  }
	}
      }
    }
    corrm_smooth[batch*num_image_elements+idx] = scale*result;
  }
}

// Smooth correlation matrices border by box filter (2D)
template<class REAL> __global__ void
smooth_correlation_matrices_border_kernel( typename complext<REAL>::Type *corrm, typename complext<REAL>::Type *corrm_smooth, intd<2>::Type image_dims, unsigned int number_of_border_threads )
{
  const int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const int batch = blockIdx.y;

  const int num_image_elements = prod(image_dims);

  if( idx < number_of_border_threads ){
    
    intd<2>::Type co;
    const int half_width = kernel_width>>1;

    co.vec[1] = idx/image_dims.vec[0];
    co.vec[1] = min(co.vec[1], half_width );
    
    if( co.vec[1] == half_width ){
      int new_idx = idx-half_width*image_dims.vec[0];
      int num_skips = new_idx/half_width;
      int rows_offset = min(num_skips>>1, image_dims.vec[1]-(half_width<<1) );
      co.vec[1] += rows_offset;

      if( co.vec[1] == (half_width + image_dims.vec[1]-(half_width<<1)) ){
	new_idx -= ((image_dims.vec[1]-(half_width<<1))*(half_width<<1));
        co.vec[1] += (new_idx / image_dims.vec[0]);
	co.vec[0] = (new_idx % image_dims.vec[0]);
      }
      else{
	co.vec[0] = (num_skips%2)*(image_dims.vec[0]-half_width) + (new_idx%half_width);
      }
    }
    else{
      co.vec[0] = idx%image_dims.vec[0];
    }
    
    const int x = co.vec[0];
    const int y = co.vec[1];
    
    const int size_x = image_dims.vec[0];
    const int size_y = image_dims.vec[1];
    
    const int yminus = y-half_width;
    const int xminus = x-half_width;

    const REAL scale = get_one<REAL>()/((REAL)(kernel_width*kernel_width));
    
    typename complext<REAL>::Type result = get_zero<typename complext<REAL>::Type>();
 
#pragma unroll
    for (int ky = 0; ky < kernel_width; ky++) {
#pragma unroll
      for (int kx = 0; kx < kernel_width; kx++) {
	
	if( (yminus+ky >=0) ){
	  if( yminus+ky < size_y ){
	    if( xminus+kx >= 0 ){
	      if( xminus+kx < size_x ){
		
		int source_offset = 
		  batch*num_image_elements +
		  (yminus+ky)*size_x +
		  (xminus+kx);
		
		result += corrm[source_offset];
	      }
	    }
	  }
	}
      }
    }
    corrm_smooth[batch*num_image_elements+co_to_idx<2>(co,image_dims)] = scale*result;  
  }
}

template<class REAL, unsigned int D> void
smooth_correlation_matrices( cuNDArray<typename complext<REAL>::Type> *corrm, cuNDArray<typename complext<REAL>::Type> *corrm_smooth )
{
  typename intd<D>::Type image_dims;

  for( unsigned int i=0; i<D; i++ ){
    image_dims.vec[i] = corrm->get_size(i);
  }
  
  unsigned int number_of_batches = 1;
  
  for( unsigned int i=D; i<corrm->get_number_of_dimensions(); i++ ){
    number_of_batches *= corrm->get_size(i);
  }
  
  int device; cudaGetDevice( &device );
  cudaDeviceProp deviceProp; cudaGetDeviceProperties( &deviceProp, device );

  dim3 blockDim(deviceProp.maxThreadsPerBlock);
  dim3 gridDim((unsigned int) ceil((double)prod(image_dims)/blockDim.x), number_of_batches);

  smooth_correlation_matrices_kernel<REAL><<<gridDim, blockDim>>>
    ( corrm->get_data_ptr(), corrm_smooth->get_data_ptr(), image_dims );
  
  CHECK_FOR_CUDA_ERROR();

  unsigned int number_of_border_threads = ((kernel_width>>1)<<1)*(sum(image_dims)-((kernel_width>>1)<<1));
  blockDim = dim3(128);
  gridDim = dim3((unsigned int) ceil((double)number_of_border_threads/blockDim.x), number_of_batches);
  
  smooth_correlation_matrices_border_kernel<REAL><<<gridDim, blockDim>>>
    ( corrm->get_data_ptr(), corrm_smooth->get_data_ptr(), image_dims, number_of_border_threads );

  CHECK_FOR_CUDA_ERROR();
}

extern __shared__ char shared_mem[];

// Extract CSM
template<class REAL> __global__ void
extract_csm_kernel( typename complext<REAL>::Type *corrm, typename complext<REAL>::Type *csm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int i = threadIdx.x;
  
  if( idx < num_elements ){    
    
    // Get the dominant eigenvector for each correlation matrix.
    // Copying Peter Kellman's approach we use the power method:
    //  b_k+1 = A*b_k / ||A*b_k||
    
    typename complext<REAL>::Type *data_out = (typename complext<REAL>::Type*) shared_mem;
    typename complext<REAL>::Type *tmp_v = &(((typename complext<REAL>::Type*) shared_mem)[num_batches*blockDim.x]);
    
    const unsigned int iterations = 2;
    
    for( unsigned int c=0; c<num_batches; c++){
      data_out[c*blockDim.x+i] = get_one<typename complext<REAL>::Type >();
    }
    
    for( unsigned int it=0; it<iterations; it++ ){
      
      for( unsigned int c=0; c<num_batches; c++){
	tmp_v[c*blockDim.x+i] = get_zero<typename complext<REAL>::Type >();
      }
      
      for( unsigned j=0; j<num_batches; j++){
	for( unsigned int k=0; k<num_batches; k++){
	  tmp_v[j*blockDim.x+i] += corrm[(k*num_batches+j)*num_elements+idx]*data_out[k*blockDim.x+i];
	}
      }
      
      REAL tmp = get_zero<REAL>();
      
      for (unsigned int c=0; c<num_batches; c++){
	tmp += norm_squared(tmp_v[c*blockDim.x+i]);
      }
      
      tmp = sqrt(tmp);
      tmp = reciprocal(tmp);
      
      for (unsigned int c=0; c<num_batches; c++){
	typename complext<REAL>::Type res = tmp*tmp_v[c*blockDim.x+i];
	data_out[c*blockDim.x+i] = res;
      }
    }
    
    for (unsigned int c=0; c<num_batches; c++){
      csm[c*num_elements+idx] = data_out[c*blockDim.x+i];
    }
  }
}

// Extract CSM
template<class REAL> __global__ void
extract_csm_kernel( typename complext<REAL>::Type *corrm, typename complext<REAL>::Type *csm, unsigned int num_batches, unsigned int num_elements, typename complext<REAL>::Type *tmp_v )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < num_elements ){    
    
    // Get the dominant eigenvector for each correlation matrix.
    // Copying Peter Kellman's approach we use the power method:
    //  b_k+1 = A*b_k / ||A*b_k||
    
    const unsigned int iterations = 2;

    for( unsigned int c=0; c<num_batches; c++){
      csm[c*num_elements+idx] = get_one<typename complext<REAL>::Type >();
    }
    
    for( unsigned int it=0; it<iterations; it++ ){

      for( unsigned int c=0; c<num_batches; c++){
	tmp_v[c*num_elements+idx] = get_zero<typename complext<REAL>::Type >();
      }
      
      for( unsigned j=0; j<num_batches; j++){
	for( unsigned int k=0; k<num_batches; k++){
	  typedef typename complext<REAL>::Type T;
	  tmp_v[j*num_elements+idx] += mul<T>(corrm[(k*num_batches+j)*num_elements+idx],csm[k*num_elements+idx]);
	}
      }

      REAL tmp = get_zero<REAL>();
      
      for (unsigned int c=0; c<num_batches; c++){
	tmp += norm_squared<REAL>(tmp_v[c*num_elements+idx]);
      }
      
      tmp = sqrt(tmp);
      tmp = reciprocal(tmp);
      
      for (unsigned int c=0; c<num_batches; c++){
	typename complext<REAL>::Type res = tmp*tmp_v[c*num_elements+idx];
	csm[c*num_elements+idx] = res;
      }
    }
  }
}

// Extract CSM
template<class REAL> __host__ 
boost::shared_ptr<cuNDArray<typename complext<REAL>::Type> > extract_csm(cuNDArray<typename complext<REAL>::Type> *corrm_in, unsigned int number_of_batches, unsigned int number_of_elements )
{
  vector<unsigned int> image_dims;

  for( unsigned int i=0; i<corrm_in->get_number_of_dimensions()-1; i++ ){
    image_dims.push_back(corrm_in->get_size(i));
  }
  
  // Allocate output
  cuNDArray<typename complext<REAL>::Type> *out = new cuNDArray<typename complext<REAL>::Type>; out->create(&image_dims);

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  /*  
      if( out != 0x0 )
      extract_csm_kernel<REAL><<< gridDim, blockDim, number_of_batches*blockDim.x*2*sizeof(typename complext<REAL>::Type) >>>
      ( corrm_in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
  */

  // Temporary buffer. TODO: use shared memory
  cuNDArray<typename complext<REAL>::Type> *tmp_v = new cuNDArray<typename complext<REAL>::Type>; tmp_v->create(&image_dims);

  if( out != 0x0 && tmp_v != 0x0 )
    extract_csm_kernel<REAL><<< gridDim, blockDim >>>
      ( corrm_in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements, tmp_v->get_data_ptr() );

  CHECK_FOR_CUDA_ERROR();
  
  delete tmp_v;
  return boost::shared_ptr<cuNDArray<typename complext<REAL>::Type> >(out);
}

// Set refence phase
template<class REAL> __global__ void
set_phase_reference_kernel( typename complext<REAL>::Type *csm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < num_elements ){
    REAL angle = arg<REAL>(csm[idx]); //Phase of the first coil
    REAL sin_a, cos_a; gad_sincos( angle, &sin_a, &cos_a );

    typename complext<REAL>::Type tmp;
    tmp.vec[0] = cos_a; tmp.vec[1] = sin_a;
    tmp = conj<typename complext<REAL>::Type>(tmp);

    for( unsigned int c=0; c<num_batches; c++ ){
      typename complext<REAL>::Type val = csm[c*num_elements+idx];
      typedef typename complext<REAL>::Type T;
      val = mul<T>( val, tmp );
      csm[c*num_elements+idx] = val;
    }
  }
}
  
// Set reference phase
template<class REAL> __host__ 
void set_phase_reference(cuNDArray<typename complext<REAL>::Type> *csm, unsigned int number_of_batches, unsigned int number_of_elements )
{
  dim3 blockDim(128);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  set_phase_reference_kernel<REAL><<< gridDim, blockDim >>>( csm->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

//
// Template instantiation
//

//template EXPORTGPUPMRI boost::shared_ptr< cuNDArray<typename complext<float>::Type > > estimate_b1_map<float,1>(cuNDArray<typename complext<float>::Type >*);
template EXPORTGPUPMRI boost::shared_ptr< cuNDArray<typename complext<float>::Type > > estimate_b1_map<float,2>(cuNDArray<typename complext<float>::Type >*);
//template boost::shared_ptr< cuNDArray<typename complext<float>::Type > > estimate_b1_map<float,3>(cuNDArray<typename complext<float>::Type >*);
//template boost::shared_ptr< cuNDArray<typename complext<float>::Type > > estimate_b1_map<float,4>(cuNDArray<typename complext<float>::Type >*);

//template EXPORTGPUPMRI boost::shared_ptr< cuNDArray<typename complext<double>::Type > > estimate_b1_map<double,1>(cuNDArray<typename complext<double>::Type >*);
template EXPORTGPUPMRI boost::shared_ptr< cuNDArray<typename complext<double>::Type > > estimate_b1_map<double,2>(cuNDArray<typename complext<double>::Type >*);
//template EXPORTGPUPMRI boost::shared_ptr< cuNDArray<typename complext<double>::Type > > estimate_b1_map<double,3>(cuNDArray<typename complext<double>::Type >*);
//template EXPORTGPUPMRI boost::shared_ptr< cuNDArray<typename complext<double>::Type > > estimate_b1_map<double,4>(cuNDArray<typename complext<double>::Type >*);
