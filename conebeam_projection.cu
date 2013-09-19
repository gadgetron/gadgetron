//
// This code performs 3D cone beam CT forwards and backwards projection
//

#include "conebeam_projection.h"
#include "float3x3.h"
#include "hoCuNDArray_operators.h"
#include "vector_td.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_utils.h"
#include "cuNFFT.h"
#include "check_CUDA.h"
#include "GPUTimer.h"
#include "cudaDeviceManager.h"
#include "hoNDArray_fileio.h"
#include "setup_grid.h"

#include <cuda_runtime_api.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#define PS_ORIGIN_CENTERING
#define IS_ORIGIN_CENTERING
//#define FLIP_Z_AXIS

// Read the projection/image data respectively as a texture (for input)
// - taking advantage of the cache and hardware interpolation
//

#define NORMALIZED_TC 1

static texture<float, 3, cudaReadModeElementType> 
image_tex( NORMALIZED_TC, cudaFilterModeLinear, cudaAddressModeBorder );

static texture<float, cudaTextureType2DLayered, cudaReadModeElementType> 
projections_tex( NORMALIZED_TC, cudaFilterModeLinear, cudaAddressModeBorder );

namespace Gadgetron 
{

  inline __host__ __device__
  float degrees2radians(float degree) {
    return degree * (CUDART_PI_F/180.0f);
  }
  
  //
  // Redundancy correction for short scan mode
  // - i.e. when you do not acquire a full rotation of data
  //
  // See "Optimal short scan convolution reconstruction for fanbeam CT", Dennis Parker, Med. Phys. 9(2) 1982
  // and (for the implementation) "Parker weights revisited", Wesarg et al, Med. Phys. 29(3) 2002.
  //
  
  static const float epsilon = 0.001f; // to avoid singularities
  
  static __inline__ __device__ float S( float beta )
  {
    if( beta <= -0.5f ) return 0.0f;
    else if( beta > -0.5f && beta < 0.5f ) return 0.5f*(1+sinf(CUDART_PI_F*beta));
    else /*if( beta >= 0.5f )*/ return 1.0f;
  }
  
  static __inline__ __device__ float B( float alpha, float delta )
  {
    return 2.0f*(delta-alpha)+epsilon;
  }
  
  static __inline__ __device__ float b( float alpha, float delta )
  {
    const float q = 0.1f; // with q=1 this formulae reduce to conventional Parker weights
    return q*B(alpha, delta);
  }
  
  __global__ void
  redundancy_correct_kernel( float *projections,
			     float *angles,
				 uintd3 dims, // Dimensions of the projections array
			     float delta // The half-fan angle
			     )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const unsigned int num_elements = prod(dims);
    
    if( idx < num_elements ){
	
      const float in = projections[idx];
      const uintd3 co = idx_to_co<3>( idx, dims );
      const float tan_delta = tanf(delta);
      const float alpha = -atanf((float(co[0])/float(dims[0])-0.5f)*2.0f*tan_delta);
      const float beta = degrees2radians(angles[co[2]]);
      
      float omega = 0.5f*(S(beta/b(alpha, delta)-0.5f)+
			  S((beta+2.0f*(alpha-delta)-epsilon)/b(alpha, delta)+0.5f)-
			  S((beta-CUDART_PI_F+2.0f*alpha)/b(-alpha, delta)-0.5f)-
			  S((beta-CUDART_PI_F-2.0f*delta-epsilon)/b(-alpha, delta)+0.5f));
      
      projections[idx] = in*omega;
    }      
  }

  void
  redundancy_correct( cuNDArray<float> *projections,
		      std::vector<float> angles,
		      float delta // The half-fan angle in radians
		      ) 
  {
    
    //
    // Validate the input
    //

    if( projections == 0x0 ){
      throw std::runtime_error("Error: redundancy_correct: illegal array pointer provided");
    }

    if( projections->get_number_of_dimensions() != 3 ){
      throw std::runtime_error("Error: redundancy_correct: projections array must be three-dimensional");
    }

    if( projections->get_size(2) != angles.size() ){
      throw std::runtime_error("Error: redundancy_correct: inconsistent sizes of input array/vector");
    }

    const unsigned int projection_res_x = projections->get_size(0);
    const unsigned int projection_res_y = projections->get_size(1);
    const unsigned int num_projections = projections->get_size(2);
    uintd3 dims(projection_res_x, projection_res_y, num_projections);

    // If not the case already, then try to transform the angles array into the range [0;360[
    //
    /*
    const float min_angle = *std::min_element(angles.begin(),angles.end());
    const float max_angle = *std::max_element(angles.begin(),angles.end());

    if( max_angle - min_angle > 360.0f ){
      throw std::runtime_error("Error: redundancy_correct: angles array range greater than 360 degrees currently not supported");
    }

    if( min_angle < -360.0f ){
      throw std::runtime_error("Error: redundancy_correct: unexpected range of the angles array (1)");
    }

    if( max_angle > 360.0f ){
      throw std::runtime_error("Error: redundancy_correct: unexpected range of the angles array (2)");
    }

    for( unsigned int i = 0; i<angles.size(); i++ ){
      angles[i] = std::fmod(angles[i]+360.0f, 360.0f);
      }*/

    float *angles_DevPtr;

    cudaMalloc( (void**) &angles_DevPtr, angles.size()*sizeof(float));
    cudaMemcpy( angles_DevPtr, &angles[0], angles.size()*sizeof(float), cudaMemcpyHostToDevice );
    CHECK_FOR_CUDA_ERROR();

    // Launch kernel
    //

    dim3 dimBlock, dimGrid;
    setup_grid( prod(dims), &dimBlock, &dimGrid );

    write_nd_array<float>( projections->to_host().get(), "before.real" );

    redundancy_correct_kernel<<< dimGrid, dimBlock >>>
      ( projections->get_data_ptr(), angles_DevPtr, dims, delta );

    write_nd_array<float>( projections->to_host().get(), "after.real" );

    CHECK_FOR_CUDA_ERROR();

    cudaFree(angles_DevPtr);
    CHECK_FOR_CUDA_ERROR();
  }


  //
  // Forwards projection
  //

  __global__ void
  conebeam_forwards_projection_kernel( float *projections,
				       float *angles,
				       floatd2 *offsets,
				       floatd3 is_dims_in_pixels,
				       floatd3 is_dims_in_mm,
				       intd2 ps_dims_in_pixels_int,
				       floatd2 ps_dims_in_mm,
				       int num_projections,
				       float SDD,
				       float SAD,
				       const int num_samples_per_ray )
  {


    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const int num_elements = prod(ps_dims_in_pixels_int)*num_projections;

    if( idx < num_elements){

      const intd3 co = idx_to_co<3>( idx, intd3(ps_dims_in_pixels_int[0], ps_dims_in_pixels_int[1], num_projections) );

      // Projection space dimensions and spacing
      //

      const floatd2 ps_dims_in_pixels = floatd2(ps_dims_in_pixels_int[0], ps_dims_in_pixels_int[1]);
      const floatd2 ps_spacing = ps_dims_in_mm / ps_dims_in_pixels;

      // Determine projection angle and rotation matrix
      //

      const float angle = angles[co[2]];
      const float3x3 rotation = calcRotationMatrixAroundZ(degrees2radians(angle));

      // Find start and end point for the line integral (image space)
      //

      floatd3 startPoint = floatd3(0.0f, -SAD, 0.0f);
      startPoint = mul(rotation, startPoint);

      // Projection plate indices
      //

#ifdef PS_ORIGIN_CENTERING
      const floatd2 ps_pc = floatd2(co[0], co[1]) + floatd2(0.5);
#else
      const floatd2 ps_pc = floatd2(co[0], co[1]);
#endif

      // Convert the projection plate coordinates into image space,
      // - local to the plate in metric units
      // - including half-fan and sag correction

      const floatd2 proj_coords = (ps_pc / ps_dims_in_pixels - 0.5f) * ps_dims_in_mm;// + offsets[co[2]];

      // Define the end point for the line integrals
      //

      const float ADD = SDD - SAD; // in mm.
      floatd3 endPoint = floatd3(proj_coords[0], ADD, proj_coords[1]);
      endPoint = mul(rotation, endPoint);

      // Find direction vector of the line integral
      //

      floatd3 dir = endPoint-startPoint;

      // Perform integration only inside the bounding cylinder of the image volume
      //
      /*
	const floatd2 is_dims_in_mm_xy(is_dims_in_mm[0], is_dims_in_mm[1]);
	const float radius = norm( is_dims_in_mm_xy * floatd2(0.5f) );
	const floatd3 unitDir = dir / SDD;

	startPoint = startPoint + (unitDir * (SAD-radius));
	endPoint = endPoint - (unitDir * (ADD-radius));
      */
      const floatd3 start = amin((-is_dims_in_mm-startPoint)/dir,(is_dims_in_mm-startPoint)/dir);
      const floatd3 end = amax((-is_dims_in_mm-startPoint)/dir,(is_dims_in_mm-startPoint)/dir);

      float a1 = fmax(max(start),0.0f);
      float aend = fmin(min(end),1.0f);
      //endPoint = startPoint+aend*dir;
      startPoint += a1*dir;

      const float sampling_distance = norm((aend-a1)*dir)/num_samples_per_ray;
      // Now perform conversion of the line integral start/end into voxel coordinates
      //

      startPoint /= is_dims_in_mm;
#ifdef FLIP_Z_AXIS
      startPoint[2] *= -1.0f;
#endif
      startPoint += 0.5f;
      //startPoint *= is_dims_in_pixels;
      dir /= is_dims_in_mm;
#ifdef FLIP_Z_AXIS
      dir[2] *= -1.0f;
#endif
      //dir *= is_dims_in_pixels;
      dir /= float(num_samples_per_ray); // now in step size units

      //
      // Perform line integration
      //

      float result = 0.0f;

      for ( int sampleIndex = 0; sampleIndex<num_samples_per_ray; sampleIndex++) {

#ifndef IS_ORIGIN_CENTERING
	floatd3 samplePoint = startPoint+dir*float(sampleIndex) + floatd3(0.5f)/is_dims_in_pixels;
#else
	floatd3 samplePoint = startPoint+dir*float(sampleIndex);
#endif

	// Accumulate result
	//

	result += tex3D( image_tex, samplePoint[0], samplePoint[1], samplePoint[2] );
      }

      // Normalization factor
      //

      // Output
      //

      projections[idx] = result*sampling_distance;
    }
  }

  //
  // Forwards projection of a 3D volume onto a set of (binned) projections
  //

  void
  conebeam_forwards_projection( hoCuNDArray<float> *projections,
				hoCuNDArray<float> *image,
				std::vector<float> angles,
				std::vector<floatd2> offsets,
				std::vector<unsigned int> indices,
				int projections_per_batch,
				float samples_per_pixel,
				floatd3 is_dims_in_mm,
				floatd2 ps_dims_in_mm,
				float SDD,
				float SAD,
				bool accumulate )
  {

    //
    // Validate the input
    //


    if( projections == 0x0 || image == 0x0 ){
      throw std::runtime_error("Error: conebeam_forwards_projection: illegal array pointer provided");
    }

    if( projections->get_number_of_dimensions() != 3 ){
      throw std::runtime_error("Error: conebeam_forwards_projection: projections array must be three-dimensional");
    }

    if( image->get_number_of_dimensions() != 3 ){
      throw std::runtime_error("Error: conebeam_forwards_projection: image array must be three-dimensional");
    }

    if( projections->get_size(2) != angles.size() || projections->get_size(2) != offsets.size() ) {
      throw std::runtime_error("Error: conebeam_forwards_projection: inconsistent sizes of input arrays/vectors");
    }

    int projection_res_x = projections->get_size(0);
    int projection_res_y = projections->get_size(1);

    int num_projections_in_bin = indices.size();
    int num_projections_in_all_bins = projections->get_size(2);

    int matrix_size_x = image->get_size(0);
    int matrix_size_y = image->get_size(1);
    int matrix_size_z = image->get_size(2);


    hoCuNDArray<float> *int_projections = projections;
    if (accumulate){
      int_projections = new hoCuNDArray<float>(projections->get_dimensions());
    }
    if( projections_per_batch > num_projections_in_bin )
      projections_per_batch = num_projections_in_bin;

    int num_batches = (num_projections_in_bin+projections_per_batch-1) / projections_per_batch;

    // Build texture from input image
    //
    cudaFuncSetCacheConfig(conebeam_forwards_projection_kernel, cudaFuncCachePreferL1);
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaExtent extent;
    extent.width = matrix_size_x;
    extent.height = matrix_size_y;
    extent.depth = matrix_size_z;

    cudaMemcpy3DParms cpy_params = {0};
    cpy_params.kind = cudaMemcpyHostToDevice;
    cpy_params.extent = extent;

    cudaArray *image_array;
    cudaMalloc3DArray(&image_array, &channelDesc, extent);
    CHECK_FOR_CUDA_ERROR();

    cpy_params.dstArray = image_array;
    cpy_params.srcPtr = make_cudaPitchedPtr
      ((void*)image->get_data_ptr(), extent.width*sizeof(float), extent.width, extent.height);
    cudaMemcpy3D(&cpy_params);
    CHECK_FOR_CUDA_ERROR();

    cudaBindTextureToArray(image_tex, image_array, channelDesc);
    CHECK_FOR_CUDA_ERROR();

    // Allocate the angles, offsets and projections in device memory
    // Using as many streams as batches enables concurrent copy/compute on devices with "just" one copy engine
    //


    float* projections_DevPtr;
    float* projections_DevPtr2;

    cudaStream_t mainStream;
    cudaStream_t indyStream;


    cudaStreamCreate(&mainStream);
    cudaStreamCreate(&indyStream);


    cudaMalloc( (void**) &projections_DevPtr, projection_res_x*projection_res_y*projections_per_batch*sizeof(float));

    cudaMalloc( (void**) &projections_DevPtr2, projection_res_x*projection_res_y*projections_per_batch*sizeof(float));

    std::vector<float> angles_vec;
    std::vector<floatd2> offsets_vec;

    for( int p=0; p<indices.size(); p++ ){

      int from_id = indices[p];

      if( from_id >= num_projections_in_all_bins ) {
	throw std::runtime_error("Error: conebeam_forwards_projection: illegal index in bin");
      }

      angles_vec.push_back(angles[from_id]);
      offsets_vec.push_back(offsets[from_id]);
    }

    thrust::device_vector<float> angles_devVec(angles_vec);
    thrust::device_vector<floatd2> offsets_devVec(offsets_vec);

    for (unsigned int batch=0; batch<num_batches; batch++ ){

      //
      // Iterate over the batches to copy input
      //

      int from_projection = batch * projections_per_batch;
      int to_projection = (batch+1) * projections_per_batch;

      if (to_projection > num_projections_in_bin)
	to_projection = num_projections_in_bin;

      int projections_in_batch = to_projection-from_projection;


      // Block/grid configuration
      //

      dim3 dimBlock, dimGrid;
      setup_grid( projection_res_x*projection_res_y*projections_in_batch, &dimBlock, &dimGrid );

      //
      // Launch kernel
      //

      floatd3 is_dims_in_pixels(matrix_size_x, matrix_size_y, matrix_size_z);
      intd2 ps_dims_in_pixels(projection_res_x, projection_res_y);



      float* raw_angles = thrust::raw_pointer_cast(&angles_devVec[from_projection]);
      floatd2* raw_offsets = thrust::raw_pointer_cast(&offsets_devVec[from_projection]);
      conebeam_forwards_projection_kernel<<< dimGrid, dimBlock, 0, mainStream >>>
	( projections_DevPtr, raw_angles, raw_offsets,
	  is_dims_in_pixels, is_dims_in_mm,
	  ps_dims_in_pixels, ps_dims_in_mm,
	  projections_in_batch, SDD, SAD, samples_per_pixel*float(matrix_size_x) );


      //If not initial batch, start copying the old stuff

      int p = from_projection;
      while( p<to_projection) {

	int num_sequential_projections = 1;
	while( p+num_sequential_projections < to_projection &&
	       indices[p+num_sequential_projections]==(indices[p+num_sequential_projections-1]+1) ){
	  num_sequential_projections++;
	}

	int to_id = indices[p];
	int size = projection_res_x*projection_res_y;


	cudaMemcpyAsync( int_projections->get_data_ptr()+to_id*size,
			 projections_DevPtr+(p-from_projection)*size,
			 size*num_sequential_projections*sizeof(float),
			 cudaMemcpyDeviceToHost, mainStream);
	//	CHECK_FOR_CUDA_ERROR();

	p += num_sequential_projections;
      }



      std::swap(projections_DevPtr,projections_DevPtr2);
      std::swap(mainStream, indyStream);

    }


    cudaFree(projections_DevPtr);
    cudaFree(projections_DevPtr2);
    cudaFreeArray(image_array);

    CHECK_FOR_CUDA_ERROR();

    if (accumulate){
      *projections += *int_projections;
      delete int_projections;
    }
  }

  __global__ void
  conebeam_backwards_projection_kernel( float *image,
					float *angles,
					floatd2 *offsets,
					intd3 is_dims_in_pixels_int,
					floatd3 is_dims_in_mm,
					floatd2 ps_dims_in_pixels,
					floatd2 ps_dims_in_mm,
					int num_projections_in_batch,
					float num_projections_in_bin,
					float SDD,
					float SAD,
					bool use_fbp,
					bool accumulate )
  {

    // Image voxel to backproject into (pixel coordinate and index)
    //

    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const int num_elements = prod(is_dims_in_pixels_int);

    if( idx < num_elements ){

      const intd3 co = idx_to_co<3>(idx, is_dims_in_pixels_int);

#ifdef IS_ORIGIN_CENTERING
      const floatd3 is_pc = floatd3(co[0], co[1], co[2]) + floatd3(0.5);
#else
      const floatd3 is_pc = floatd3(co[0], co[1], co[2]);
#endif

      // Normalized image space coordinate [-0.5, 0.5[
      //

      const floatd3 is_dims_in_pixels(is_dims_in_pixels_int[0],is_dims_in_pixels_int[1],is_dims_in_pixels_int[2]);

#ifdef FLIP_Z_AXIS
      floatd3 is_nc = is_pc / is_dims_in_pixels - floatd3(0.5f);
      is_nc[2] *= -1.0f;
#else
      const floatd3 is_nc = is_pc / is_dims_in_pixels - floatd3(0.5f);
#endif

      // Image space coordinate in metric units
      //

      const floatd3 pos = is_nc * is_dims_in_mm;

      // Read the existing output value for accumulation at this point.
      // The cost of this fetch is hidden by the loop

      const float incoming = (accumulate) ? image[idx] : 0.0f;

      // Backprojection loop
      //

      float result = 0.0f;

      for( int projection = 0; projection < num_projections_in_batch; projection++ ) {

	// Projection angle
	const float angle = degrees2radians(angles[projection]);

	// Projection rotation matrix
	const float3x3 inverseRotation = calcRotationMatrixAroundZ(-angle);

	// Rotated image coordinate (local to the projection's coordinate system)
	const floatd3 pos_proj = mul(inverseRotation, pos);

	// Project the image position onto the projection plate.
	// Account for half-fan and sag offsets.

	const floatd3 startPoint = floatd3(0.0f, -SAD, 0.0f);
	floatd3 dir = pos_proj - startPoint;
	dir = dir / dir[1];
	const floatd3 endPoint = startPoint + dir * SDD;
	const floatd2 endPoint2d = floatd2(endPoint[0], endPoint[2]);// - offsets[projection];

	// Convert metric projection coordinates into pixel coordinates
	//

#ifndef PS_ORIGIN_CENTERING
	floatd2 ps_pc = ((endPoint2d / ps_dims_in_mm) + floatd2(0.5f)) + floatd2(0.5f)/ps_dims_in_pixels;
	//floatd2 ps_pc = ((endPoint2d / ps_dims_in_mm) + floatd2(0.5f)) * ps_dims_in_pixels + floatd2(0.5f);
#else
	floatd2 ps_pc = ((endPoint2d / ps_dims_in_mm) + floatd2(0.5f));
#endif

	// Apply filter (filtered backprojection mode only)
	//

	float weight = 1.0;

	if( use_fbp ){

	  // Equation 3.59, page 96 and equation 10.2, page 386
	  // in Computed Tomography 2nd edition, Jiang Hsieh
	  //

	  const float xx = pos[0];
	  const float yy = pos[1];
	  const float beta = angle;
	  const float r = ::hypot(xx,yy);
	  const float phi = ::atan2(yy,xx);
	  const float D = SAD;
	  const float ym = r*::sin(beta-phi);
	  const float U = (D+ym)/D;
	  weight = 1.0f/(U*U);
	}

	// Read the projection data (bilinear interpolation enabled) and accumulate
	//


	result +=  weight * tex2DLayered( projections_tex, ps_pc[0], ps_pc[1], projection );

      }

      // Output normalized image
      //

      image[idx] = incoming + result / num_projections_in_bin;
    }
  }

  // Compute cosine weights (for filtered backprojection)
  //

  boost::shared_ptr< cuNDArray<float> >
  get_cosinus_weights( intd2 ps_dims_in_pixels, floatd2 ps_dims_in_mm, float D0d, float Ds0 )
  {
    std::vector<unsigned int> dims;
    dims.push_back(ps_dims_in_pixels[0]);
    dims.push_back(ps_dims_in_pixels[1]);

    hoCuNDArray<float> weights(&dims);
    float* data = weights.get_data_ptr();

    const float Dsd = D0d + Ds0;

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( int y=0; y<ps_dims_in_pixels[1]; y++ ) {
      for( int x=0; x<ps_dims_in_pixels[0]; x++ ) {

	double xx = (( double(x) / double(ps_dims_in_pixels[0])) - 0.5) * ps_dims_in_mm[0];
	double yy = (( double(y) / double(ps_dims_in_pixels[1])) - 0.5) * ps_dims_in_mm[1];
	double s = Ds0 * xx/Dsd;
	double t = Ds0 * yy/Dsd;

	// Equation 10.1, page 386 in Computed Tomography 2nd edition, Jiang Hsieh
	//

	double value = Ds0 / std::sqrt( Ds0*Ds0 + s*s + t*t );
	data[x+y*ps_dims_in_pixels[0]] = float(value);
      }
    }
    return boost::shared_ptr< cuNDArray<float> >(new cuNDArray<float>(&weights));
  }

  //
  // Use ramp filter in filtered backprojection
  //

  boost::shared_ptr< cuNDArray<float> >
  get_ramp( int dim, float delta )
  {
    std::vector<unsigned int> dims, dims_os;
    dims.push_back(dim);
    dims_os.push_back(dim<<1);

    hoCuNDArray<float> weights(&dims);
    float* data = weights.get_data_ptr();

    // Compute ramp weights
    //

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( int i=0; i<dim; i++ ) {

      // Equation 3.29, page 73 in Computed Tomography 2nd edition, Jiang Hsieh
      //

      int n = i-dim/2;
      float value;
      if(n==0) {
	value = 1.0f/(4.0f*delta*delta);
      } else if ((i%2)==0) { // even
	value = 0.0f;
      } else { // odd
	float tmp = n*CUDART_PI_F*delta;
	value = -1.0f/(tmp*tmp);
      }
      data[i] = value;
    }

    cuNDArray<float> tmp_real(&weights);
    cuNDArray<float> weights_os(&dims_os);
    pad<float,1>( &tmp_real, &weights_os );

    // Compute FFT of the ramp weights
    //

    boost::shared_ptr< cuNDArray<float_complext> > tmp_cplx = real_to_complex<float_complext>(&weights_os);
    cuNFFT_plan<float,1>().fft(tmp_cplx.get(), cuNFFT_plan<float,1>::NFFT_FORWARDS);
    boost::shared_ptr< cuNDArray<float> >  tmp_res = real(tmp_cplx.get());

    cuNDArray<float> *res = new cuNDArray<float>(&dims);
    crop<float,1>( from_std_vector<unsigned int,1>(dims)>>1, tmp_res.get(), res );

    return boost::shared_ptr< cuNDArray<float> >(res);
  }

  //
  // Backprojection
  //

  void
  conebeam_backwards_projection( hoCuNDArray<float> *projections,
				 hoCuNDArray<float> *image,
				 std::vector<float> angles,
				 std::vector<floatd2> offsets,
				 std::vector<unsigned int> indices,
				 int projections_per_batch,
				 intd3 is_dims_in_pixels,
				 floatd3 is_dims_in_mm,
				 floatd2 ps_dims_in_mm,
				 float SDD,
				 float SAD,
				 bool use_fbp,
				 bool use_oversampling_in_fbp,
				 float maximum_angle,
				 bool accumulate )
  {

    //
    // Validate the input
    //

    if( projections == 0x0 || image == 0x0 ){
      throw std::runtime_error("Error: conebeam_backwards_projection: illegal array pointer provided");
    }

    if( projections->get_number_of_dimensions() != 3 ){
      throw std::runtime_error("Error: conebeam_backwards_projection: projections array must be three-dimensional");
    }

    if( image->get_number_of_dimensions() != 3 ){
      throw std::runtime_error("Error: conebeam_backwards_projection: image array must be three-dimensional");
    }

    if( projections->get_size(2) != angles.size() || projections->get_size(2) != offsets.size() ) {
      throw std::runtime_error("Error: conebeam_backwards_projection: inconsistent sizes of input arrays/vectors");
    }

    // Some utility variables
    //

    int matrix_size_x = image->get_size(0);
    int matrix_size_y = image->get_size(1);
    int matrix_size_z = image->get_size(2);

    unsigned int warp_size = cudaDeviceManager::Instance()->warp_size();

    if( matrix_size_x < warp_size ){
      throw std::runtime_error("Error: conebeam_backwards_projection: image dimension 'x' must be at least as large as the warp size");
    }

    if( matrix_size_x % warp_size ){
      throw std::runtime_error("Error: conebeam_backwards_projection: image dimension 'x' must be a multiplum of the warp size");
    }

    floatd3 is_dims(matrix_size_x, matrix_size_y, matrix_size_z);
    int num_image_elements = matrix_size_x*matrix_size_y*matrix_size_z;

    int projection_res_x = projections->get_size(0);
    int projection_res_y = projections->get_size(1);

    floatd2 ps_dims_in_pixels(projection_res_x, projection_res_y);

    int num_projections_in_all_bins = projections->get_size(2);
    int num_projections_in_bin = indices.size();

    if( projections_per_batch > num_projections_in_bin )
      projections_per_batch = num_projections_in_bin;

    int num_batches = (num_projections_in_bin+projections_per_batch-1) / projections_per_batch;

    // Allocate device memory for the backprojection result and constant memory arrays
    //

    boost::shared_ptr< cuNDArray<float> > image_device;

    if( accumulate ){
      image_device = boost::shared_ptr< cuNDArray<float> >
	(new cuNDArray<float>(image));
    }
    else{
      image_device = boost::shared_ptr< cuNDArray<float> >
	(new cuNDArray<float>(image->get_dimensions().get()));
    }

    float *angles_DevPtr;
    Gadgetron::floatd2 *offsets_DevPtr;

    cudaMalloc( (void**) &angles_DevPtr, projections_per_batch*sizeof(float) );
    CHECK_FOR_CUDA_ERROR();

    cudaMalloc( (void**) &offsets_DevPtr, projections_per_batch*sizeof(floatd2) );
    CHECK_FOR_CUDA_ERROR();

    //
    // Iterate over batches
    //

    for( int batch = 0; batch < num_batches; batch++ ) {

      int from_projection = batch * projections_per_batch;
      int to_projection = (batch+1) * projections_per_batch;

      if (to_projection > num_projections_in_bin )
	to_projection = num_projections_in_bin;

      int projections_in_batch = to_projection-from_projection;

      
	printf("batch: %03i, handling projections: %03i - %03i, angles: %.2f - %.2f\n", 
	batch, from_projection, to_projection-1, angles[from_projection], angles[to_projection-1]);      
      

      // Allocate device memory for projections and upload
      //
      std::vector<unsigned int> dims;
      dims.push_back(projection_res_x);
      dims.push_back(projection_res_y);
      dims.push_back(projections_in_batch);

      cuNDArray<float> projections_batch(&dims);

      // Copy result from device to host adhering to the binning
      // Be sure to copy sequentially numbered projections in one copy operation.

      int p = from_projection;
      while( p<to_projection ) {

	int num_sequential_projections = 1;
	while( p+num_sequential_projections < to_projection &&
	       indices[p+num_sequential_projections]==(indices[p+num_sequential_projections-1]+1) ){
	  num_sequential_projections++;
	}

	int from_id = indices[p];
	int size = projection_res_x*projection_res_y;

	cudaMemcpy( projections_batch.get_data_ptr()+(p-from_projection)*size,
		    projections->get_data_ptr()+from_id*size,
		    size*num_sequential_projections*sizeof(float), cudaMemcpyHostToDevice );

	CHECK_FOR_CUDA_ERROR();

	p += num_sequential_projections;
      }

      // Make vectors of the angles and offsets present in the current bin (indices)
      //

      std::vector<float> angles_vec;
      std::vector<floatd2> offsets_vec;

      for( int p=from_projection; p<to_projection; p++ ){

	int from_id = indices[p];

	if( from_id >= num_projections_in_all_bins ) {
	  throw std::runtime_error("Error: conebeam_backwards_projection: illegal index in bin");
	}

	angles_vec.push_back(angles[from_id]);
	offsets_vec.push_back(offsets[from_id]);
      }

      // Upload angles and offsets data to device
      //

      cudaMemcpy(angles_DevPtr, &angles_vec[0], projections_in_batch*sizeof(float), cudaMemcpyHostToDevice);
      CHECK_FOR_CUDA_ERROR();

      cudaMemcpy(offsets_DevPtr, &offsets_vec[0], projections_in_batch*sizeof(floatd2), cudaMemcpyHostToDevice);
      CHECK_FOR_CUDA_ERROR();

      // Pre-backprojection filtering (if enabled)
      //

      if (use_fbp) {

	static boost::shared_ptr< cuNDArray<float> > cos_weights;
	static bool initialized = false;

	// Compute cosine weighing at first invocation
	//

	if (!initialized) {
	  const float ADD = SDD - SAD;
	  cos_weights = get_cosinus_weights(intd2(projection_res_x, projection_res_y), ps_dims_in_mm, ADD, SAD);
	  initialized = true;
	  //write_nd_array<float>( cos_weights->to_host().get(), "cos_weights.real");
	}

	// 1. Cosine weighting "SDD / sqrt(SDD*SDD + u*u + v*v)"
	//

	projections_batch *= *cos_weights;

	// 1.5 redundancy correct (in short scan mode)
	//

	if( maximum_angle < 360.0f ){
	  float delta = std::atan(ps_dims_in_mm[0]/(2.0f*SDD));
	  redundancy_correct( &projections_batch, angles_vec, delta );
	}

	// 2. Apply ramp filter (using FFT)
	//

	// Use oversamping if asked for...
	//

	if( use_oversampling_in_fbp ) {

	  std::vector<unsigned int> dims_os;
	  dims_os.push_back(projection_res_x<<1);
	  dims_os.push_back(projection_res_y);
	  dims_os.push_back(projections_in_batch);

	  cuNDArray<float> projections_os(&dims_os);
	  pad<float,2>( &projections_batch, &projections_os );

	  static boost::shared_ptr< cuNDArray<float> > ramp;
	  static bool ramp_initialized = false;

	  if (!ramp_initialized) {
	    ramp = get_ramp( projection_res_x<<1, ps_dims_in_mm[0]/(projection_res_x<<1) );
	    ramp_initialized = true;
	    //write_nd_array<float>( ramp->to_host().get(), "ramp.real");
	  }

	  boost::shared_ptr< cuNDArray<complext<float> > > complex_projections_os =
	    real_to_complex< complext<float> >(&projections_os);
	  cuNFFT_plan<float,1>().fft(complex_projections_os.get(),cuNFFT_plan<float,1>::NFFT_FORWARDS);
	  *complex_projections_os *= *ramp;
	  cuNFFT_plan<float,1>().fft(complex_projections_os.get(),cuNFFT_plan<float,1>::NFFT_BACKWARDS);
	  projections_os = *real(complex_projections_os.get());

	  uintd<2>::Type offset;
	  offset.vec[0]= projection_res_x>>1;
	  offset.vec[1]= 0;
	  crop<float,2>( offset, &projections_os, &projections_batch );
	}
	else{

	  static boost::shared_ptr< cuNDArray<float> > ramp;
	  static bool ramp_initialized = false;

	  if (!ramp_initialized) {
	    ramp = get_ramp( projection_res_x, ps_dims_in_mm[0]/projection_res_x );
	    ramp_initialized = true;
	    //write_nd_array<float>( ramp->to_host().get(), "ramp.real");
	  }

	  boost::shared_ptr< cuNDArray<complext<float> > > complex_projections =
	    real_to_complex< complext<float> >(&projections_batch);
	  cuNFFT_plan<float,1>().fft(complex_projections.get(),cuNFFT_plan<float,1>::NFFT_FORWARDS);
	  *complex_projections *= *ramp;
	  cuNFFT_plan<float,1>().fft(complex_projections.get(),cuNFFT_plan<float,1>::NFFT_BACKWARDS);
	  projections_batch = *real(complex_projections.get());
	}
      }

      // Build array for input texture
      //

      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      cudaExtent extent;
      extent.width = projection_res_x;
      extent.height = projection_res_y;
      extent.depth = projections_in_batch;

      cudaArray *projections_array;
      cudaMalloc3DArray( &projections_array, &channelDesc, extent, cudaArrayLayered );
      CHECK_FOR_CUDA_ERROR();

      cudaMemcpy3DParms cpy_params = {0};
      cpy_params.extent = extent;
      cpy_params.dstArray = projections_array;
      cpy_params.kind = cudaMemcpyDeviceToDevice;
      cpy_params.srcPtr =
	make_cudaPitchedPtr( (void*)projections_batch.get_data_ptr(), projection_res_x*sizeof(float),
			     projection_res_x, projection_res_y );
      cudaMemcpy3D( &cpy_params );
      CHECK_FOR_CUDA_ERROR();

      cudaBindTextureToArray( projections_tex, projections_array, channelDesc );
      CHECK_FOR_CUDA_ERROR();

      // Define dimensions of grid/blocks.
      dim3 dimBlock, dimGrid;
      setup_grid( matrix_size_x*matrix_size_y*matrix_size_z, &dimBlock, &dimGrid );

      //
      // Go...
      //

      cudaFuncSetCacheConfig(conebeam_backwards_projection_kernel, cudaFuncCachePreferL1);

      // Invoke kernel
      conebeam_backwards_projection_kernel<<< dimGrid, dimBlock >>>
	( image_device->get_data_ptr(), angles_DevPtr, offsets_DevPtr,
	  is_dims_in_pixels, is_dims_in_mm, ps_dims_in_pixels, ps_dims_in_mm,
	  projections_in_batch, num_projections_in_bin,
	  SDD, SAD, use_fbp,
	  (batch==0) ? accumulate : true );

      CHECK_FOR_CUDA_ERROR();

      /*
	static int counter = 0;
	char filename[256];
	sprintf((char*)filename, "bp_%d.real", counter);
	write_nd_array<float>( image_device->to_host().get(), filename );
	counter++; */

      // Cleanup
      cudaUnbindTexture(projections_tex);
      cudaFreeArray(projections_array);
      CHECK_FOR_CUDA_ERROR();
    }

    // Copy result from device to host
    cudaMemcpy( image->get_data_ptr(), image_device->get_data_ptr(),
		num_image_elements*sizeof(float), cudaMemcpyDeviceToHost );
    CHECK_FOR_CUDA_ERROR();

    cudaFree(offsets_DevPtr);
    cudaFree(angles_DevPtr);
    CHECK_FOR_CUDA_ERROR();
  }
}
