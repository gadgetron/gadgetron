//
// This code performs 3D cone beam CT forwards and backwards projection
//

#include "conebeam_projection.h"
#include "float3x3.h"
#include "hoNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_utils.h"
#include "cuNFFT.h"
#include "check_CUDA.h"
#include "GPUTimer.h"
#include "GadgetronException.h"
#include "cudaDeviceManager.h"
#include "hoNDArray_fileio.h"

#include <cuda_runtime_api.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#define PS_ORIGIN_CENTERING
#define IS_ORIGIN_CENTERING
//#define FLIP_Z_AXIS

// Cuda insists on that constant variables and textures are defined at the global scope.
// This WILL become problematic in a multi-threaded (multi-GPU) environment one day.
// Acces should really be protected by mutexes.

static /*__constant__*/ /*__device__*/ float *angles_DevPtr;
static /*__constant__*/ /*__device__*/ Gadgetron::floatd2 *offsets_DevPtr;

#define NORMALIZED_TC 0

static texture<float, 3, cudaReadModeElementType> 
image_tex( NORMALIZED_TC, cudaFilterModeLinear, cudaAddressModeBorder );

static texture<float, cudaTextureType2DLayered, cudaReadModeElementType> 
projections_tex( NORMALIZED_TC, cudaFilterModeLinear, cudaAddressModeBorder );

namespace Gadgetron 
{
  
  // Read the projection/image data respectively as a texture (taking advantage of the cache and interpolation)
  //

  inline __host__ __device__
  float degrees2radians(float degree) {
    return degree * (CUDART_PI_F/180.0f);
  }

  //
  // Redundancy correction is used for FDK in short fan mode (when you do not acquire a full rotation of data)
  //

  __global__ void
  redundancy_correct_kernel( float *projection,
			     float *angles,
			     float lambda_max, 
			     float gamma_m, 
			     float gamma_thresh )
  {
    const int PROJ_IDX_X = threadIdx.x+blockIdx.x*blockDim.x;
    const int PROJ_IDX_Y = blockIdx.y;
    const int PROJ_ID = blockIdx.z;
    const int PROJECTION_RES_X = blockDim.x*gridDim.x;
    const int PROJECTION_RES_Y = gridDim.y;

    const long int idx = PROJ_ID*PROJECTION_RES_Y*PROJECTION_RES_X*PROJ_IDX_Y*PROJECTION_RES_X+PROJ_IDX_X;
    const float gamma = -gamma_m + 2.0f*gamma_m*(float(PROJ_IDX_X)/float(PROJECTION_RES_X));
    const float lambda = angles[PROJ_ID]*(CUDART_PI_F/180.0f);

    float res = 0.0f;

    if( lambda >= 0.0f && lambda < 2.0f*(gamma_thresh+gamma) ){
      float tmp = ::sin(CUDART_PI_F/4.0f*lambda/(gamma_thresh+gamma));
      res = tmp*tmp;
    }
    else if( lambda >= 2.0f*(gamma_thresh+gamma) && lambda < CUDART_PI_F+2.0f*gamma ){
      res = 1.0f;
    }
    else if( lambda >= CUDART_PI_F+2.0f*gamma && lambda < CUDART_PI_F+2.0f*gamma_thresh ){
      float tmp = ::sin(CUDART_PI_F/4.0f*(CUDART_PI_F+2.0f*gamma_thresh-lambda)/(gamma_thresh+gamma));
      res = tmp*tmp;
    }
    else ;

    projection[idx] *= res;
  }

  void
  redundancy_correct( cuNDArray<float> *projections,
		      std::vector<float> angles,
		      float lambda_max, 
		      float gamma_m, 
		      float gamma_thresh )
  {
    
    //
    // Validate the input 
    //
    
    if( projections == 0x0 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: redundancy_correct: illegal array pointer provided"));
    }
    
    if( projections->get_number_of_dimensions() != 3 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: redundancy_correct: projections array must be three-dimensional"));
    }
    
    if( projections->get_size(2) != angles.size() ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: redundancy_correct: inconsistent sizes of input array/vector"));
    }

    const int projection_res_x = projections->get_size(0);
    const int projection_res_y = projections->get_size(1);
    const int num_projections = projections->get_size(2);

    // If not the case already, then try to transform the angles array into the range [0;360[
    //

    const float min_angle = *std::min_element(angles.begin(),angles.end());

    if( min_angle < -180.0f ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: redundancy_correct: unexpected range of the angles array"));
    }
    
    if( min_angle < 0.0f ){
      transform(angles.begin(), angles.end(), angles.begin(), bind2nd(std::plus<double>(), 180.0f ));
    }

    const float max_angle = *std::max_element(angles.begin(),angles.end());
    if( max_angle*(CUDART_PI_F/180.0f) >= CUDART_PI_F+2.0f*gamma_thresh ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: redundancy_correct: unexpected range of the angles array after correction"));
    }

    cudaMalloc( (void**) &angles_DevPtr, angles.size()*sizeof(float));
    cudaMemcpy( angles_DevPtr, &angles[0], angles.size()*sizeof(float), cudaMemcpyHostToDevice );
    CHECK_FOR_CUDA_ERROR();

    // We have been somewhat lazy in the reundancy kernel (and the block/grid setup). 
    // If the checks below fail, it is a fairly easy fix...
    //

    int warp_size = cudaDeviceManager::Instance()->warp_size();
  
    if( projection_res_x < warp_size ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: redundancy_correct: projection sizes have to be at least as large as the warp size"));
    }

    if( projection_res_x % warp_size ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: redundancy_correct: projection sizes have to be a multiplum of the warp size"));
    }

    // Launch kernel
    //

    dim3 dimBlock(std::min(projection_res_x, 256));
    dim3 dimGrid((projection_res_x+dimBlock.x-1)/dimBlock.x, projection_res_y, num_projections);
    
    redundancy_correct_kernel<<< dimGrid, dimBlock >>>
      ( projections->get_data_ptr(), angles_DevPtr, lambda_max, gamma_m, gamma_thresh );

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
				       floatd3 is_spacing_in_mm,
				       intd2 ps_dims_in_pixels_int, 
				       floatd2 ps_dims_in_mm,
				       float SDD, 
				       float SAD,
				       float num_samples_per_ray, 
				       bool accumulate )
  {
    // Some defines to give the thread/block/grid setup (more) meaningful names
    //

    const int PROJ_IDX_X = threadIdx.x + blockIdx.x * blockDim.x;
    const int PROJ_IDX_Y = threadIdx.y + blockIdx.y * blockDim.y;
    const int PROJECTION = blockIdx.z;
    const int PROJECTION_RES_X = ps_dims_in_pixels_int[0];
    const int PROJECTION_RES_Y = ps_dims_in_pixels_int[1];

    if ( PROJ_IDX_X < PROJECTION_RES_X && PROJ_IDX_Y < PROJECTION_RES_Y) {

      // Projection space pixel index
      //

      const long int projectionSpaceId =
	PROJ_IDX_X +
	PROJ_IDX_Y * PROJECTION_RES_X +
	PROJECTION * PROJECTION_RES_X * PROJECTION_RES_Y;

      // Projection space dimensions and spacing
      //

      const floatd2 ps_dims_in_pixels = floatd2(PROJECTION_RES_X, PROJECTION_RES_Y);
      const floatd2 ps_spacing = ps_dims_in_mm / ps_dims_in_pixels;

      // Determine projection angle and rotation matrix
      //

      const float angle = angles[PROJECTION];
      const float3x3 rotation = calcRotationMatrixAroundZ(degrees2radians(angle));

      // Find start and end point for the line integral (image space)
      //

      floatd3 startPoint = floatd3(0.0f, -SAD, 0.0f);
      startPoint = mul(rotation, startPoint);

      // Projection plate indices
      //

#ifdef PS_ORIGIN_CENTERING
      const floatd2 ps_pc = floatd2(PROJ_IDX_X, PROJ_IDX_Y) + floatd2(0.5);
#else
      const floatd2 ps_pc = floatd2(PROJ_IDX_X, PROJ_IDX_Y);
#endif

      // Convert the projection plate coordinates into image space,
      // - local to the plate in metric units
      // - including half-fan and sag correction

      const floatd2 proj_coords = (ps_pc / ps_dims_in_pixels - 0.5f) * ps_dims_in_mm + offsets[PROJECTION];

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

      const floatd3 is_dims_in_mm = is_dims_in_pixels*is_spacing_in_mm;
      const floatd2 is_dims_in_mm_xy(is_dims_in_mm[0], is_dims_in_mm[1]);
      const floatd2 is_spacing_in_mm_xy(is_spacing_in_mm[0], is_spacing_in_mm[1]);
      const float radius = norm( (is_dims_in_mm_xy + is_spacing_in_mm_xy) * floatd2(0.5f) );
      const floatd3 unitDir = dir / SDD;

      startPoint = startPoint + (unitDir * (SAD-radius));
      endPoint = endPoint - (unitDir * (ADD-radius));

      // Now perform conversion of the line integral start/end into voxel coordinates
      //      

      startPoint /= is_dims_in_mm;
#ifdef FLIP_Z_AXIS
      startPoint[2] *= -1.0f;
#endif
      startPoint += 0.5f;
      startPoint *= is_dims_in_pixels;
      dir /= is_dims_in_mm;
#ifdef FLIP_Z_AXIS
      dir[2] *= -1.0f;
#endif
      dir *= is_dims_in_pixels;
      dir /= num_samples_per_ray; // now in step size units

      // Read the existing output value for accumulation at this point.
      // The cost of this fetch is hidden from our subsequent integration
    
      const float incoming = (accumulate) ? projections[projectionSpaceId] : 0.0f;

      // Normalization factor
      // 

      const float rayLength = norm(startPoint-endPoint)/num_samples_per_ray;

      //
      // Perform line integration
      //

      float result = 0.0f;

      for ( float sampleIndex = 0.0f; sampleIndex<num_samples_per_ray-0.5f; sampleIndex+= 1.0f) {
      
#ifndef IS_ORIGIN_CENTERING
	floatd3 samplePoint = startPoint+dir*sampleIndex + floatd3(0.5f);
#else
	floatd3 samplePoint = startPoint+dir*sampleIndex;
#endif

	// Accumulate result
	//

	result += tex3D( image_tex, samplePoint[0], samplePoint[1], samplePoint[2] );
      }

      // And output
      //

      projections[projectionSpaceId] = incoming + result*rayLength;
    }
  }

  //
  // Forwards projection of a 3D volume onto a set of (binned) projections
  //

  void 
  conebeam_forwards_projection( hoNDArray<float> *projections,
				hoNDArray<float> *image,
				std::vector<float> angles, 
				std::vector<floatd2> offsets, 
				std::vector<unsigned int> indices,
				int projections_per_batch, 
				int num_samples_per_ray,
				floatd3 is_spacing_in_mm, 
				floatd2 ps_dims_in_mm,
				float SDD, 
				float SAD,
				bool accumulate )
  {

    //
    // Validate the input 
    //

    if( projections == 0x0 || image == 0x0 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_forwards_projection: illegal array pointer provided"));
    }

    if( projections->get_number_of_dimensions() != 3 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_forwards_projection: projections array must be three-dimensional"));
    }

    if( image->get_number_of_dimensions() != 3 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_forwards_projection: image array must be three-dimensional"));
    }

    if( projections->get_size(2) != angles.size() || projections->get_size(2) != offsets.size() ) {
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_forwards_projection: inconsistent sizes of input arrays/vectors"));
    }

    int projection_res_x = projections->get_size(0);
    int projection_res_y = projections->get_size(1);

    int num_projections_in_bin = indices.size();
    int num_projections_in_all_bins = projections->get_size(2);

    int matrix_size_x = image->get_size(0);
    int matrix_size_y = image->get_size(1);
    int matrix_size_z = image->get_size(2);

    if( projections_per_batch > num_projections_in_bin )
      projections_per_batch = num_projections_in_bin;

    int num_batches = (num_projections_in_bin+projections_per_batch-1) / projections_per_batch;

    // Build texture from input image
    //

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

    // Keep the angles and offsets lists in constant device memory
    //
  
    cudaMalloc( (void**) &angles_DevPtr, projections_per_batch*sizeof(float));
    CHECK_FOR_CUDA_ERROR();
  
    cudaMalloc( (void**) &offsets_DevPtr, projections_per_batch*sizeof(floatd2));
    CHECK_FOR_CUDA_ERROR();
  
    // Allocate device memory for the projections
    //

    float* projections_DevPtr;
    cudaMalloc( (void**) &projections_DevPtr, projection_res_x*projection_res_y*projections_per_batch*sizeof(float));
    CHECK_FOR_CUDA_ERROR();

    //
    // Iterate over the batches
    //
  
    for( int batch = 0; batch < num_batches; batch++ ) {

      int from_projection = batch * projections_per_batch;
      int to_projection = (batch+1) * projections_per_batch;

      if (to_projection > num_projections_in_bin)
	to_projection = num_projections_in_bin;

      int projections_in_batch = to_projection-from_projection;

      printf("batch: %03i, handling projections: %03i - %03i, angle %.2f - %.2f\n",
	     batch, from_projection, to_projection-1,
	     angles[from_projection], angles[to_projection-1]);

      // Make vectors of the angles and offsets present in the current bin (indices)
      //
      
      std::vector<float> angles_vec;
      std::vector<floatd2> offsets_vec;

      for( int p=from_projection; p<to_projection; p++ ){

	int from_id = indices[p];

	if( from_id >= num_projections_in_all_bins ) {
	  BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_forwards_projection: illegal index in bin"));
	}
	
	angles_vec.push_back(angles[from_id]);
	offsets_vec.push_back(offsets[from_id]);
      }
	
      // Upload angles and offsets data to device
      //
      
      cudaMemcpy(angles_DevPtr, &angles_vec[0], projections_in_batch*sizeof(float), cudaMemcpyHostToDevice);
      CHECK_FOR_CUDA_ERROR();
    
      cudaMemcpy(offsets_DevPtr, &offsets_vec[0], projections_in_batch*sizeof(float), cudaMemcpyHostToDevice);
      CHECK_FOR_CUDA_ERROR();

      // Block/grid configuration
      //

      dim3 dimBlock(std::min(projection_res_x, 256));
      dim3 dimGrid((projection_res_x+dimBlock.x-1)/dimBlock.x, 
		   projection_res_y, 
		   projections_in_batch );

      //
      // Launch kernel 
      //

      floatd3 is_dims_in_pixels(matrix_size_x, matrix_size_y, matrix_size_z);
      intd2 ps_dims_in_pixels(projection_res_x, projection_res_y);
      
      cudaFuncSetCacheConfig(conebeam_forwards_projection_kernel, cudaFuncCachePreferL1);

      conebeam_forwards_projection_kernel<<< dimGrid, dimBlock >>>
	( projections_DevPtr, angles_DevPtr, offsets_DevPtr,
	  is_dims_in_pixels, is_spacing_in_mm, 
	  ps_dims_in_pixels, ps_dims_in_mm,
	  SDD, SAD,
	  num_samples_per_ray, 
	  false );
      
      CHECK_FOR_CUDA_ERROR();

      // Copy result from device to host adhering to the binning
      // Be sure to copy sequentially numbered projections in one copy operation.
      
      int p = from_projection;
      while( p<to_projection ) {
	
	int num_sequential_projections = 1;
	while( p+num_sequential_projections < to_projection && 
	       indices[p+num_sequential_projections]==(indices[p+num_sequential_projections-1]+1) ){
	  num_sequential_projections++;
	}
	
	int to_id = indices[p];
	int size = projection_res_x*projection_res_y;

	if( !accumulate ) {
	  cudaMemcpy( projections->get_data_ptr()+to_id*size, 
		      projections_DevPtr+(p-from_projection)*size, 
		      size*num_sequential_projections*sizeof(float), cudaMemcpyDeviceToHost );
	}
	else {
	  std::vector<unsigned int> dims;
	  dims.push_back(projection_res_x);
	  dims.push_back(projection_res_y);
	  dims.push_back(num_sequential_projections);
	  hoNDArray<float> tmp(&dims);
	  cudaMemcpy( tmp.get_data_ptr(), projections_DevPtr+(p-from_projection)*size, 
		      size*num_sequential_projections*sizeof(float), cudaMemcpyDeviceToHost );
	  hoNDArray<float> target( &dims, projections->get_data_ptr()+to_id*size );
	  target += tmp;
	}
	
	CHECK_FOR_CUDA_ERROR();
	
	p += num_sequential_projections;
      }
    }

    // Clean up
    //

    cudaFree(projections_DevPtr);
    cudaFree(offsets_DevPtr);
    cudaFree(angles_DevPtr);
    cudaUnbindTexture(image_tex);
    cudaFreeArray(image_array);

    CHECK_FOR_CUDA_ERROR();
  }

  __global__ void 
  conebeam_backwards_projection_kernel( float *image, 
					float *angles,
					floatd2 *offsets,
					floatd3 is_dims_in_pixels, 
					floatd3 is_spacing_in_mm,
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

    const int IMAGE_X = threadIdx.x;
    const int IMAGE_Y = blockIdx.x;
    const int IMAGE_Z = blockIdx.y;

    const long int id =
      IMAGE_X +
      IMAGE_Y * blockDim.x +
      IMAGE_Z * blockDim.x*gridDim.x;

#ifdef IS_ORIGIN_CENTERING
    const floatd3 is_pc = floatd3( IMAGE_X, IMAGE_Y, IMAGE_Z ) + floatd3(0.5);
#else
    const floatd3 is_pc = floatd3( IMAGE_X, IMAGE_Y, IMAGE_Z );
#endif

    // Normalized image space coordinate [-0.5, 0.5[
    //

#ifdef FLIP_Z_AXIS
    floatd3 is_nc = is_pc / is_dims_in_pixels - floatd3(0.5f);
    is_nc[2] *= -1.0f;
#else
    const floatd3 is_nc = is_pc / is_dims_in_pixels - floatd3(0.5f);
#endif
    
    // Image space coordinate in metric units
    //

    const floatd3 is_dims_in_mm = is_spacing_in_mm * is_dims_in_pixels;
    const floatd3 pos = is_nc * is_dims_in_mm;

    // Read the existing output value for accumulation at this point.
    // The cost of this fetch is hidden by the loop
    
    const float incoming = (accumulate) ? image[id] : 0.0f;
    
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
      const floatd2 endPoint2d = floatd2(endPoint[0], endPoint[2]) - offsets[projection];

      // Convert metric projection coordinates into pixel coordinates
      //

#ifndef PS_ORIGIN_CENTERING
      floatd2 ps_pc = ((endPoint2d / ps_dims_in_mm) + floatd2(0.5f)) * ps_dims_in_pixels + floatd2(0.5f);
#else
      floatd2 ps_pc = ((endPoint2d / ps_dims_in_mm) + floatd2(0.5f)) * ps_dims_in_pixels;
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

      const float val = weight * tex2DLayered( projections_tex, ps_pc[0], ps_pc[1], projection );
      result += val;
    }
    image[id] = incoming + result / num_projections_in_bin;
  }

  // Compute cosine weights (for filtered backprojection)
  //

  boost::shared_ptr< cuNDArray<float> > 
  get_cosinus_weights( intd2 ps_dims_in_pixels, floatd2 ps_dims_in_mm, float D0d, float Ds0 ) 
  {
    std::vector<unsigned int> dims;
    dims.push_back(ps_dims_in_pixels[0]);
    dims.push_back(ps_dims_in_pixels[1]);

    hoNDArray<float> weights(&dims);
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

    hoNDArray<float> weights(&dims);
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
  conebeam_backwards_projection( hoNDArray<float> *projections, 
				 hoNDArray<float> *image,
				 std::vector<float> angles, 
				 std::vector<floatd2> offsets, 
				 std::vector<unsigned int> indices,
				 int projections_per_batch,
				 uintd3 is_dims_in_pixels, 
				 floatd3 is_spacing_in_mm, 
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
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_backwards_projection: illegal array pointer provided"));
    }

    if( projections->get_number_of_dimensions() != 3 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_backwards_projection: projections array must be three-dimensional"));
    }

    if( image->get_number_of_dimensions() != 3 ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_backwards_projection: image array must be three-dimensional"));
    }

    if( projections->get_size(2) != angles.size() || projections->get_size(2) != offsets.size() ) {
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_backwards_projection: inconsistent sizes of input arrays/vectors"));
    }

    // Some utility variables
    //
    
    int matrix_size_x = image->get_size(0);
    int matrix_size_y = image->get_size(1);
    int matrix_size_z = image->get_size(2);

    unsigned int warp_size = cudaDeviceManager::Instance()->warp_size();
  
    if( matrix_size_x < warp_size ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_backwards_projection: image dimension 'x' must be at least as large as the warp size"));
    }

    if( matrix_size_x % warp_size ){
      BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_backwards_projection: image dimension 'x' must be a multiplum of the warp size"));
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

      printf("batch: %03i, handling projections: %03i - %03i, angle %.2f - %.2f\n", batch,
	     from_projection, to_projection-1, angles[from_projection], angles[to_projection-1]);      

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
	  BOOST_THROW_EXCEPTION(runtime_error("Error: conebeam_backwards_projection: illegal index in bin"));
	}
	
	angles_vec.push_back(angles[from_id]);
	offsets_vec.push_back(offsets[from_id]);	
      }
      
      // Upload angles and offsets data to device
      //
      
      cudaMemcpy(angles_DevPtr, &angles_vec[0], projections_in_batch*sizeof(float), cudaMemcpyHostToDevice);
      CHECK_FOR_CUDA_ERROR();
    
      cudaMemcpy(offsets_DevPtr, &offsets_vec[0], projections_in_batch*sizeof(float), cudaMemcpyHostToDevice);
      CHECK_FOR_CUDA_ERROR();

      // Pre-backprojection filtering (if enabled)
      //

      if (use_fbp) {	

	static boost::shared_ptr< cuNDArray<float> > cos_weights; 
	static bool initialized = false;

	// Compute cosine weighing at first invocation
	if (!initialized) {
	  const float ADD = SDD - SAD;
	  cos_weights = get_cosinus_weights(intd2(projection_res_x, projection_res_y), ps_dims_in_mm, ADD, SAD);
	  initialized = true;
	  //write_nd_array<float>( cos_weights->to_host().get(), "cos_weights.real");
	}

	// 1. Cosine weighting "SDD / sqrt( SDD*SDD + u*u + v*v)"
	//   

	projections_batch *= *cos_weights;

	// 1.5 redundancy correct (for short fan scans)
	//

	if( maximum_angle < 2.0*CUDART_PI_F ){
	  float gamma_m = CUDART_PI_F/2.0-::atan(SDD/(ps_dims_in_mm[0]/2.0f));
	  float gamma_thresh = (maximum_angle-CUDART_PI_F)/2.0f;
	  redundancy_correct( &projections_batch, angles_vec, maximum_angle, gamma_m, gamma_thresh );
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
      dim3 dimBlock( matrix_size_x );
      dim3 dimGrid( matrix_size_y, matrix_size_z );

      //
      // Go...
      //
      
      cudaFuncSetCacheConfig(conebeam_backwards_projection_kernel, cudaFuncCachePreferL1);

      // Invoke kernel
      conebeam_backwards_projection_kernel<<< dimGrid, dimBlock >>>
	( image_device->get_data_ptr(), angles_DevPtr, offsets_DevPtr,
	  floatd3(float(is_dims_in_pixels[0]), float(is_dims_in_pixels[1]), float(is_dims_in_pixels[2])),
	  is_spacing_in_mm, ps_dims_in_pixels, ps_dims_in_mm,
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
