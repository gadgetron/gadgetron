#include "b1_map.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "vector_td_utilities.h"
#include "real_utilities.h"
#include "real_utilities_device.h"
#include "complext.h"
#include "check_CUDA.h"
#include "cudaDeviceManager.h"
#include "setup_grid.h"

#include <iostream>
#include <cmath>

using namespace std;

namespace Gadgetron{

  const int kernel_width = 7;

  template<class REAL, unsigned int D> static void smooth_correlation_matrices( cuNDArray<complext<REAL> >*, cuNDArray<complext<REAL> >*);
  template<class REAL> static cuNDArray<complext<REAL>> extract_csm(const cuNDArray<complext<REAL> >&, unsigned int, unsigned int);
  template<class REAL> static void set_phase_reference( cuNDArray<complext<REAL> >*, unsigned int, unsigned int);
  template<class T> static void find_stride( cuNDArray<T> *in, unsigned int dim, unsigned int *stride, std::vector<size_t> *dims );
  template<class T> static boost::shared_ptr< cuNDArray<T> > correlation( cuNDArray<T> *in );
  template<class T> static void rss_normalize( cuNDArray<T> *in_out, unsigned int dim );

  //
  // Main method
  //

  template<class REAL, unsigned int D> cuNDArray<complext<REAL> >
  estimate_b1_map( const cuNDArray<complext<REAL> >& data_in, int target_coils)
  {

    if( data_in.get_number_of_dimensions() < 2 ){
      throw std::runtime_error("estimate_b1_map:: dimensionality mismatch.");
    }

    if( data_in.get_number_of_dimensions()-1 != D ){
      throw std::runtime_error("estimate_b1_map:: dimensionality mismatch.");
    }

    int target_coils_int = 0;
    if ((target_coils <= 0) || (target_coils > data_in.get_size(D))) {
      target_coils_int = data_in.get_size(D);
    } else {
      target_coils_int = target_coils;
    }

    vector<unsigned int> image_dims, dims_to_xform;
    unsigned int pixels_per_coil = 1;

    for( unsigned int i=0; i<D; i++ ){
      image_dims.push_back(data_in.get_size(i));
      dims_to_xform.push_back(i);
      pixels_per_coil *= data_in.get_size(i);
    }

    unsigned int ncoils = data_in.get_size(D);

    // Make a copy of input data, but only the target coils
      std::vector<size_t> odims = data_in.get_dimensions();
      odims[D] = target_coils_int;
      auto data_out = cuNDArray<complext<REAL> >(odims);

      //Now copy one coil at a time
      unsigned int elements_per_coil = data_in.get_number_of_elements()/ncoils;
      for (unsigned int i = 0; i < target_coils_int; i++) {
	cudaMemcpy(data_out.get_data_ptr()+i*elements_per_coil,
		   data_in.get_data_ptr()+i*elements_per_coil,
		   elements_per_coil*sizeof(complext<REAL>),
		   cudaMemcpyDeviceToDevice);
      }
      ncoils = target_coils_int;

    // Normalize by the RSS of the coils
    rss_normalize( &data_out, D );

    // Now calculate the correlation matrices
    boost::shared_ptr<cuNDArray<complext<REAL> > > corrm = correlation( &data_out );
    data_out.clear();

    // Smooth (onto copy of corrm)
    std::vector<size_t> dim;
    corrm->get_dimensions(dim);
    auto corrm_smooth = boost::make_shared<cuNDArray<complext<REAL>>>(dim);

    smooth_correlation_matrices<REAL,D>( corrm.get(), corrm_smooth.get() );
    corrm.reset();

    // Get the dominant eigenvector for each correlation matrix.
    auto csm = extract_csm<REAL>( *corrm_smooth, ncoils, pixels_per_coil );
    corrm_smooth.reset();

    // Set phase according to reference (coil 0)
    set_phase_reference<REAL>( &csm, ncoils, pixels_per_coil );

    return csm;
  }

  template<class T> static void find_stride( cuNDArray<T> *in, unsigned int dim,
					     unsigned int *stride, std::vector<size_t> *dims )
  {
    *stride = 1;
    for( unsigned int i=0; i<in->get_number_of_dimensions(); i++ ){
      if( i != dim )
	dims->push_back(in->get_size(i));
      if( i < dim )
	*stride *= in->get_size(i);
    }
  }

  template<class REAL, class T> __inline__  __device__ static REAL
  _rss( unsigned int idx, T *in, unsigned int stride, unsigned int number_of_batches )
  {
    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
    REAL rss = REAL(0);

    for( unsigned int i=0; i<number_of_batches; i++ )
      rss += norm(in[i*stride+in_idx]);

    rss = std::sqrt(rss);

    return rss;
  }

  template<class T> __global__ static void
  rss_normalize_kernel( T *in_out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
  {
    typedef typename realType<T>::Type REAL;

    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

    if( idx < number_of_elements ){

      REAL reciprocal_rss = 1/(_rss<REAL,T>(idx, in_out, stride, number_of_batches));

      unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);

      for( unsigned int i=0; i<number_of_batches; i++ ) {
	T out = in_out[i*stride+in_idx];
	out *= reciprocal_rss; // complex-scalar multiplication (element-wise operator)
	in_out[i*stride+in_idx] = out;
      }
    }
  }

  // Normalized RSS
  template<class T> static
  void rss_normalize( cuNDArray<T> *in_out, unsigned int dim )
  {
    unsigned int number_of_batches = in_out->get_size(dim);
    unsigned int number_of_elements = in_out->get_number_of_elements()/number_of_batches;

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( number_of_elements, &blockDim, &gridDim );

    // Find element stride
    unsigned int stride; std::vector<size_t> dims;
    find_stride<T>( in_out, dim, &stride, &dims );

    // Invoke kernel
    rss_normalize_kernel<T><<< gridDim, blockDim >>>( in_out->get_data_ptr(), stride, number_of_batches, number_of_elements );

    CHECK_FOR_CUDA_ERROR();
  }

  template<class REAL, class T> __global__ static void
  correlation_kernel( const T * __restrict__ in, T * __restrict__ corrm, unsigned int num_batches, unsigned int num_elements )
  {
    const unsigned int p = blockIdx.x*blockDim.x + threadIdx.x;
    const unsigned int i = threadIdx.y;

    if( p < num_elements ){
      for( unsigned int j=0; j<i; j++){
	T tmp = in[i*num_elements+p]*conj(in[j*num_elements+p]);
	corrm[(j*num_batches+i)*num_elements+p] = tmp;
	corrm[(i*num_batches+j)*num_elements+p] = conj(tmp);
      }
      T tmp = in[i*num_elements+p];
      corrm[(i*num_batches+i)*num_elements+p] = tmp*conj(tmp);
    }
  }

  // Build correlation matrix
  template<class T> static boost::shared_ptr< cuNDArray<T> > correlation( cuNDArray<T> *in )
  {
    typedef typename realType<T>::Type REAL;
    // Prepare internal array
    int cur_device = cudaDeviceManager::Instance()->getCurrentDevice();

    unsigned int number_of_batches = in->get_size(in->get_number_of_dimensions()-1);
    unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

    int warp_size = cudaDeviceManager::Instance()->warp_size(cur_device);
    int max_blockdim = cudaDeviceManager::Instance()->max_blockdim(cur_device);
    dim3 blockDim(((max_blockdim/number_of_batches)/warp_size)*warp_size, number_of_batches);

    if( blockDim.x == 0 ){
      throw std::runtime_error("correlation: correlation dimension exceeds device capacity.");
    }

    dim3 gridDim((number_of_elements+blockDim.x-1)/blockDim.x);

    // Invoke kernel
    std::vector<size_t> dims = in->get_dimensions(); dims.push_back(number_of_batches);
    boost::shared_ptr< cuNDArray<T> > out( new cuNDArray<T> );
    out->create(dims);

    correlation_kernel<REAL,T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );

    CHECK_FOR_CUDA_ERROR();

    return out;
  }

  // Smooth correlation matrices by box filter (1D)
  template<class REAL> __global__ static void
  smooth_correlation_matrices_kernel( const complext<REAL> * __restrict__ corrm, complext<REAL> * __restrict__ corrm_smooth, intd<1>::Type image_dims )
  {
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    const int batch = blockIdx.y;

    const int num_image_elements = prod(image_dims);

    if( idx < num_image_elements ){

      const int co = idx;
      const int x = co;

      const int size_x = image_dims.vec[0];

      const REAL scale = REAL(1)/((REAL)kernel_width);

      complext<REAL> result = complext<REAL>(0);

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
  template<class REAL> __global__ static  void
  smooth_correlation_matrices_kernel( const complext<REAL> * __restrict__ corrm, complext<REAL> * __restrict__ corrm_smooth, intd<2>::Type image_dims )
  {
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    const int batch = blockIdx.y;

    const int num_image_elements = prod(image_dims);

    if( idx < num_image_elements ){

      const intd2 co = idx_to_co(idx, image_dims);

      const int x = co.vec[0];
      const int y = co.vec[1];

      const int size_x = image_dims.vec[0];
      const int size_y = image_dims.vec[1];

      const int half_width = kernel_width>>1;

      const int yminus = y-half_width;
      const int xminus = x-half_width;
      const int yplus = y+half_width;
      const int xplus = x+half_width;

      const REAL scale = REAL(1)/((REAL)(kernel_width*kernel_width));

      complext<REAL> result = complext<REAL>(0);

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
  template<class REAL> __global__ static  void
  smooth_correlation_matrices_kernel( const  complext<REAL> * __restrict__ corrm, complext<REAL> * __restrict__ corrm_smooth, intd<3>::Type image_dims )
  {
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    const int batch = blockIdx.y;

    const int num_image_elements = prod(image_dims);

    if( idx < num_image_elements ){

      const intd3 co = idx_to_co(idx, image_dims);

      const int x = co.vec[0];
      const int y = co.vec[1];
      const int z = co.vec[2];

      const int size_x = image_dims.vec[0];
      const int size_y = image_dims.vec[1];
      const int size_z = image_dims.vec[2];

      const REAL scale = REAL(1)/((REAL)(kernel_width*kernel_width*kernel_width));

      complext<REAL> result = complext<REAL>(0);

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
  template<class REAL> __global__ static void
  smooth_correlation_matrices_kernel( const complext<REAL> * __restrict__ corrm, complext<REAL> * __restrict__ corrm_smooth, intd<4>::Type image_dims )
  {
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    const int batch = blockIdx.y;

    const int num_image_elements = prod(image_dims);

    if( idx < num_image_elements ){

      const intd4 co = idx_to_co(idx, image_dims);

      const int x = co.vec[0];
      const int y = co.vec[1];
      const int z = co.vec[2];
      const int w = co.vec[3];

      const int size_x = image_dims.vec[0];
      const int size_y = image_dims.vec[1];
      const int size_z = image_dims.vec[2];
      const int size_w = image_dims.vec[3];

      const REAL scale = REAL(1)/((REAL)(kernel_width*kernel_width*kernel_width*kernel_width));

      complext<REAL> result = complext<REAL>(0);

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

  __device__ int _min( int A, int B ){
    return (A<B) ? A : B;
  }

  // Smooth correlation matrices border by box filter (2D)
  template<class REAL> __global__ static void
  smooth_correlation_matrices_border_kernel( const complext<REAL> * __restrict__ corrm, complext<REAL> * __restrict__ corrm_smooth, intd<2>::Type image_dims, unsigned int number_of_border_threads )
  {
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    const int batch = blockIdx.y;

    const int num_image_elements = prod(image_dims);

    if( idx < number_of_border_threads ){

      intd2 co;
      const int half_width = kernel_width>>1;

      co.vec[1] = idx/image_dims.vec[0];
      co.vec[1] = _min(co.vec[1], half_width );

      if( co.vec[1] == half_width ){
        int new_idx = idx-half_width*image_dims.vec[0];
        int num_skips = new_idx/half_width;
        int rows_offset = _min(num_skips>>1, image_dims.vec[1]-(half_width<<1) );
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

      const REAL scale = REAL(1)/((REAL)(kernel_width*kernel_width));

      complext<REAL> result = complext<REAL>(0);

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
      corrm_smooth[batch*num_image_elements+co_to_idx(co,image_dims)] = scale*result;
    }
  }

  // Smooth correlation matrices border by box filter (3D)
  template<class REAL> __global__ static void
  smooth_correlation_matrices_border_kernel( const complext<REAL> * __restrict__ corrm, complext<REAL> * __restrict__ corrm_smooth, intd<3>::Type image_dims, unsigned int number_of_border_threads )
  {
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    const int batch = blockIdx.y;

    const int half_width = kernel_width >> 1;

    const int num_image_elements = prod(image_dims);
    const int num_slice_elements = image_dims.vec[0] * image_dims.vec[1];
    const int num_slice_border_elements = num_slice_elements - (image_dims.vec[0] - (half_width << 1)) * (image_dims.vec[1] - (half_width << 1));

    if (idx < number_of_border_threads) {

      intd3 co;

      co.vec[2] = idx / num_slice_elements;
      co.vec[2] = _min(co.vec[2], half_width);

      if (co.vec[2] == half_width) {
        int new_idx = idx - half_width * num_slice_elements;
        int skips = new_idx / (num_slice_border_elements >> 1);
        int offset = _min(skips >> 1, image_dims.vec[2] - (half_width << 1));
        co.vec[2] += offset;

        if (co.vec[2] == image_dims.vec[2] - half_width) {
          new_idx -= (image_dims.vec[2] - (half_width << 1)) * num_slice_border_elements;
          co.vec[2] += new_idx / num_slice_elements;
          co.vec[1] = (new_idx / image_dims.vec[0]) % image_dims.vec[1];
          co.vec[0] = new_idx % image_dims.vec[0];
        }
        else {
          new_idx %= num_slice_border_elements;
          co.vec[1] = new_idx / image_dims.vec[0];
          co.vec[1] = _min(co.vec[1], half_width);

          if (co.vec[1] == half_width) {
            new_idx -= half_width * image_dims.vec[0];
            skips = new_idx / half_width;
            offset = _min(skips >> 1, image_dims.vec[1] - (half_width << 1));
            co.vec[1] += offset;

            if (co.vec[1] == image_dims.vec[1] - half_width) {
              new_idx -= ((image_dims.vec[1] - (half_width << 1)) * (half_width << 1));
              co.vec[1] += new_idx / image_dims.vec[0];
              co.vec[0] = new_idx % image_dims.vec[0];
            }
            else {
              co.vec[0] = (skips % 2) * (image_dims.vec[0] - half_width) + (new_idx % half_width);
            }
          }
          else {
            co.vec[0] = new_idx % image_dims.vec[0];
          }
        }
      }
      else {
        co.vec[1] = (idx / image_dims.vec[0]) % image_dims.vec[1];
        co.vec[0] = idx % image_dims.vec[0];
      }

      const int x = co.vec[0];
      const int y = co.vec[1];
      const int z = co.vec[2];

      const int size_x = image_dims.vec[0];
      const int size_y = image_dims.vec[1];
      const int size_z = image_dims.vec[2];

      const int zminus = z-half_width;
      const int yminus = y-half_width;
      const int xminus = x-half_width;

      const REAL scale = REAL(1)/((REAL)(kernel_width*kernel_width*kernel_width));

      complext<REAL> result = complext<REAL>(0);

#pragma unroll
      for (int kz = 0; kz < kernel_width; kz++) {
#pragma unroll
        for (int ky = 0; ky < kernel_width; ky++) {
#pragma unroll
	        for (int kx = 0; kx < kernel_width; kx++) {
            if (zminus+kz >= 0) {
              if (zminus+kz < size_z) {
                if (yminus+ky >= 0) {
                  if (yminus+ky < size_y) {
                    if (xminus+kx >= 0) {
                      if (xminus+kx < size_x) {

                        int source_offset =
                          batch*num_image_elements +
                          (zminus+kz)*size_x*size_y +
                          (yminus+ky)*size_x +
                          (xminus+kx);

                        result += corrm[source_offset];
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      corrm_smooth[batch*num_image_elements+co_to_idx(co,image_dims)] = scale*result;
    }
  }

  template<class REAL, unsigned int D> static void
  smooth_correlation_matrices( cuNDArray<complext<REAL> > * corrm, cuNDArray<complext<REAL> > * corrm_smooth )
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
    dim3 gridDim((unsigned int) std::ceil((double)prod(image_dims)/blockDim.x), number_of_batches);

    smooth_correlation_matrices_kernel<REAL><<<gridDim, blockDim>>>
      ( corrm->get_data_ptr(), corrm_smooth->get_data_ptr(), image_dims );

    CHECK_FOR_CUDA_ERROR();
    unsigned int number_of_border_threads = prod(image_dims) - prod(image_dims-((kernel_width>>1)<<1));
    blockDim = dim3(128);
    gridDim = dim3((unsigned int) std::ceil((double)number_of_border_threads/blockDim.x), number_of_batches);

    smooth_correlation_matrices_border_kernel<REAL><<<gridDim, blockDim>>>
      ( corrm->get_data_ptr(), corrm_smooth->get_data_ptr(), image_dims, number_of_border_threads );

    CHECK_FOR_CUDA_ERROR();
  }

  extern __shared__ char shared_mem[];

  // Extract CSM
  template<class REAL> __global__ static void
  extract_csm_kernel( const complext<REAL> * __restrict__ corrm, complext<REAL> * __restrict__ csm, unsigned int num_batches, unsigned int num_elements )
  {
    const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
    const unsigned int i = threadIdx.x;

    if( idx < num_elements ){

      // Get the dominant eigenvector for each correlation matrix.
      // Copying Peter Kellman's approach we use the power method:
      //  b_k+1 = A*b_k / ||A*b_k||

      complext<REAL> *data_out = (complext<REAL>*) shared_mem;
      complext<REAL> *tmp_v = &(((complext<REAL>*) shared_mem)[num_batches*blockDim.x]);

      const unsigned int iterations = 2;

      for( unsigned int c=0; c<num_batches; c++){
	data_out[c*blockDim.x+i] = complext<REAL>(1);
      }

      for( unsigned int it=0; it<iterations; it++ ){

	for( unsigned int c=0; c<num_batches; c++){
	  tmp_v[c*blockDim.x+i] = complext<REAL>(0);
	}

	for( unsigned j=0; j<num_batches; j++){
	  for( unsigned int k=0; k<num_batches; k++){
	    tmp_v[j*blockDim.x+i] += corrm[(k*num_batches+j)*num_elements+idx]*data_out[k*blockDim.x+i];
	  }
	}

	REAL tmp = REAL(0);

	for (unsigned int c=0; c<num_batches; c++){
	  tmp += norm(tmp_v[c*blockDim.x+i]);
	}

	tmp = 1/std::sqrt(tmp);


	for (unsigned int c=0; c<num_batches; c++){
	  complext<REAL> res = tmp*tmp_v[c*blockDim.x+i];
	  data_out[c*blockDim.x+i] = res;
	}
      }

      for (unsigned int c=0; c<num_batches; c++){
	csm[c*num_elements+idx] = data_out[c*blockDim.x+i];
      }
    }
  }

  // Extract CSM
  template<class REAL> __global__ static void
  extract_csm_kernel( const complext<REAL> * __restrict__ corrm, complext<REAL> * __restrict__ csm, unsigned int num_batches, unsigned int num_elements, complext<REAL> * __restrict__ tmp_v )
  {
    const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if( idx < num_elements ){

      // Get the dominant eigenvector for each correlation matrix.
      // Copying Peter Kellman's approach we use the power method:
      //  b_k+1 = A*b_k / ||A*b_k||

      const unsigned int iterations = 2;

      for( unsigned int c=0; c<num_batches; c++){
	csm[c*num_elements+idx] = complext<REAL>(1);
      }

      for( unsigned int it=0; it<iterations; it++ ){

	for( unsigned int c=0; c<num_batches; c++){
	  tmp_v[c*num_elements+idx] = complext<REAL>(0);
	}

	for( unsigned j=0; j<num_batches; j++){
	  for( unsigned int k=0; k<num_batches; k++){
	    typedef complext<REAL> T;
	    tmp_v[j*num_elements+idx] += corrm[(k*num_batches+j)*num_elements+idx]*csm[k*num_elements+idx];
	  }
	}

	REAL tmp = REAL(0);

	for (unsigned int c=0; c<num_batches; c++){
	  tmp += norm(tmp_v[c*num_elements+idx]);
	}

	tmp = 1/std::sqrt(tmp);


	for (unsigned int c=0; c<num_batches; c++){
	  complext<REAL> res = tmp*tmp_v[c*num_elements+idx];
	  csm[c*num_elements+idx] = res;
	}
      }
    }
  }

  // Extract CSM
  template<class REAL> __host__ static
  cuNDArray<complext<REAL>> extract_csm(const cuNDArray<complext<REAL> >& corrm_in, unsigned int number_of_batches, unsigned int number_of_elements )
  {
    vector<size_t> image_dims;

    for( unsigned int i=0; i<corrm_in.get_number_of_dimensions()-1; i++ ){
      image_dims.push_back(corrm_in.get_size(i));
    }

    // Allocate output
    cuNDArray<complext<REAL> > out = cuNDArray<complext<REAL> >(image_dims);

    dim3 blockDim(256);
    dim3 gridDim((unsigned int) std::ceil((double)number_of_elements/blockDim.x));

    cuNDArray<complext<REAL> > tmp_v = cuNDArray<complext<REAL> >(image_dims);

      extract_csm_kernel<REAL><<< gridDim, blockDim >>>
	( corrm_in.get_data_ptr(), out.get_data_ptr(), number_of_batches, number_of_elements, tmp_v.get_data_ptr() );

    CHECK_FOR_CUDA_ERROR();

    return out;
  }

  // Set refence phase
  template<class REAL> __global__ static void
  set_phase_reference_kernel( complext<REAL> *csm, unsigned int num_batches, unsigned int num_elements )
  {
    const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if( idx < num_elements ){
      REAL angle = arg<REAL>(csm[idx]); //Phase of the first coil
      REAL sin_a, cos_a; gad_sincos( angle, &sin_a, &cos_a );

      complext<REAL> tmp;
      tmp._real = cos_a; tmp._imag = sin_a;
      tmp = conj(tmp);

      for( unsigned int c=0; c<num_batches; c++ ){
	complext<REAL> val = csm[c*num_elements+idx];
	typedef complext<REAL> T;
	val = val*tmp;
	csm[c*num_elements+idx] = val;
      }
    }
  }

  // Set reference phase
  template<class REAL> __host__ static
  void set_phase_reference(cuNDArray<complext<REAL> > *csm, unsigned int number_of_batches, unsigned int number_of_elements )
  {
    dim3 blockDim(128);
    dim3 gridDim((unsigned int) std::ceil((double)number_of_elements/blockDim.x));

    set_phase_reference_kernel<REAL><<< gridDim, blockDim >>>( csm->get_data_ptr(), number_of_batches, number_of_elements );

    CHECK_FOR_CUDA_ERROR();
  }



  //
  // Template instantiation
  //

  //template boost::shared_ptr< cuNDArray<complext<float> > > estimate_b1_map<float,1>(cuNDArray<complext<float> >*, int);
  template cuNDArray<complext<float>> estimate_b1_map<float,2>(const cuNDArray<complext<float> >&, int);
  template cuNDArray<complext<float>> estimate_b1_map<float,3>(const cuNDArray<complext<float> >&, int);
  //template boost::shared_ptr< cuNDArray<complext<float> > > estimate_b1_map<float,4>(cuNDArray<complext<float> >*, int);

  //template boost::shared_ptr< cuNDArray<complext<double> > > estimate_b1_map<double,1>(cuNDArray<complext<double> >*, int);
  template cuNDArray<complext<double>> estimate_b1_map<double,2>(const cuNDArray<complext<double>>&, int);
  template cuNDArray<complext<double>> estimate_b1_map<double,3>(const cuNDArray<complext<double>>&, int);
  //template boost::shared_ptr< cuNDArray<complext<double> > > estimate_b1_map<double,4>(cuNDArray<complext<double> >*, int);
}
