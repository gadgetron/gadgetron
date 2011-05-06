#include "preprocess_sense.hcu"
#include "preprocess.hcu"
#include "preprocess_private.hcu"
// TODO: How to avoid this "radial" include?
#include "preprocess_radial.hcu"
#include "FFT.hcu"
#include "NFFT.hcu"
#include "NSense_private.hcu"
#include "image_utilities.hcu"
#include "uint_util.hcu"
#include "float_util.hcu"

#include <cuComplex.h>
#include <cublas.h>
#include <math_constants.h>
#include <vector_functions.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <stdio.h>

// Utility

__global__ void 
__divide_images_kernel( unsigned int num_elements, cuFloatComplex *target, cuFloatComplex *_source, unsigned int num_targets  )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){

		cuFloatComplex source = _source[idx];
		if( cuCrealf(source) > 1.0f || cuCimagf(source) > 1.0f )
			for( unsigned int i=0; i<num_targets; i++ )
				target[i*num_elements+idx] = cuCdivf( target[i*num_elements+idx], source );
	}
}

__host__ void 
__divide_images( unsigned int num_elements, cuFloatComplex *targetDevPtr, cuFloatComplex *sourceDevPtr, unsigned int num_target_images )
{
	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/deviceProp.maxThreadsPerBlock), 1, 1 );

	// Invoke kernel
	__divide_images_kernel<<< dimGrid, dimBlock >>>( num_elements, targetDevPtr, sourceDevPtr, num_target_images );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in '__divide_images_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


/*
	Noise scaling.
	Output a dataset matching the data samples size of noise magnitudes squared outside shutter and 0 inside.
*/

// NOTICE: There is a separat version of 'noise_scaling_kernel' in preprocess_radial.cu
//		   ! Any changes made here is likely necessary in the radial case also !

__global__ void
noise_scaling_kernel( unsigned int num_points, unsigned int num_coils, cuFloatComplex *data, float4 *co, float *out_magnitudes_squared, float *out_count, float2 image_center, float shutter_radius_squared )
{
	const unsigned int sample_no = (blockIdx.x*blockDim.x+threadIdx.x);

	if( sample_no >= num_points )
		return;

	// Get sample position. 
	// We assume 2D kt-Sense for now. Time is third dimension and 'co' is stored as float4.

	float4 tmp = co[sample_no];
	float2 sample_pos = make_float2( tmp.x, tmp.y );
	float2 offset = sample_pos-image_center;

	float dist_sq = dot(offset, offset);

	// Write to global memory

	if( dist_sq >= shutter_radius_squared ){
		for( unsigned int c=0; c<num_coils; c++ ){
			cuFloatComplex _data = data[c*num_points+sample_no];
			out_magnitudes_squared[c*num_points+sample_no] = cuCrealf(cuCmulf(_data, cuConjf(_data)));
		}
		out_count[sample_no] = 1.0f;
	}
	else{
		for( unsigned int c=0; c<num_coils; c++ )
			out_magnitudes_squared[c*num_points+sample_no] = 0.0f;
		out_count[sample_no] = 0.0f;
	}
}

bool 
noise_scaling( unsigned int num_points, unsigned int num_coils, cuFloatComplex *dataDevPtr, float4 *coDevPtr, float2 image_center, float shutter_radius_squared )
{
	printf("\nERROR: Check code in 'noise_scaling'. Needs update?\n");
	return false;
/*
	if( coDevPtr == 0x0 ){
		printf("\nnoise_scaling: Trajectory device pointer is 0x0!\n");
		return false;
	}

	// First find noise magnitudes
	
	dim3 blockDim(512,1,1);
	dim3 gridDim((unsigned int) ceil((double)(num_points/blockDim.x)), 1, 1 );
	assert(num_points>=blockDim.x);
	float *noise_modulus_sq, *point_count;
	CUDA_SAFE_CALL( cudaMalloc( (void**) &noise_modulus_sq, num_coils*num_points*sizeof(float) ));
	CUDA_SAFE_CALL( cudaMalloc( (void**) &point_count, num_points*sizeof(float) ));

	noise_scaling_kernel<<< gridDim, blockDim >>>( num_points, num_coils, dataDevPtr, coDevPtr, noise_modulus_sq, point_count, image_center, shutter_radius_squared );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'noise_scaling_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Then sum for each coil and normalize

	float *noise_variances = (float*) malloc(num_coils*sizeof(float));
	float used_points = cublasSasum( num_points, point_count, 1 );

	if( used_points > 0.1f ){
		for( unsigned int c=0; c<num_coils; c++ ){
			noise_variances[c] = cublasSasum( num_points, &noise_modulus_sq[c*num_points], 1 );	
			noise_variances[c] /= used_points;
		}
	}
	else{
		printf("\nnoise_scaling: No points found within k-space shutter!\n");
		return false;
	}

	// And finally scale the samples values accordingly

	for( unsigned int c=0; c<num_coils; c++){
//		printf("\nNoise scaling with a factor of %f for coil %d", 1.0f/sqrtf(noise_variances[c]), c );
		cublasSscal( 2*num_points, 1.0f/sqrtf(noise_variances[c]), (float*)(&dataDevPtr[c*num_points]), 1 );
	}

	free(noise_variances);
	CUDA_SAFE_CALL( cudaFree( noise_modulus_sq )); noise_modulus_sq = 0x0;
	CUDA_SAFE_CALL( cudaFree( point_count )); point_count = 0x0;
	
	return true;

*/
}

template< char TYPE >
float compute_training_data_shutter_radius( unsigned int projections_per_frame )
{
/*
	if( TYPE == 1 ){

		// Find largest Fibonacci number smaller than 'number_of_projections' (F2)
		unsigned int F = 2, F1 = 1, F2 = 1, k = 1;
		float delta, gamma = (sqrtf(5.0f) + 1.0f)/2.0f, g = 1.0f/gamma;
		;

		while( F < number_of_projections ){

			F1 = F2;
			F2 = F;
			F = F1 + F2;
			k++;
		}

		unsigned int i = k>>1;

		if( k%2 == 0)
			delta = g*powf( 1.0f-g, i-1 )*M_PI;
		else
			delta = powf( 1.0f-g, i )*M_PI;

		return 1.0f / delta;
	}
	else{
		printf("\nImplementation error: 'compute_training_data_shutter_radius'. Quitting.\n");
		exit(1);
	}
*/	
	// For normal radial the formula is simply: shutter_theta = (float)projections_per_frame/M_PI;
	
	// The code above seems to be too conservative...
	return (float)projections_per_frame/(float)M_PI;
}

// Upload CSM from host
template< class UINTd, class FLOATd, char TYPE > bool
upload_csm( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE> *plan, cuFloatComplex *csmHostPtr )
{
	/*
		Some initial checks for input/plan consistency
	*/

	if( !plan->initialized || (plan->number_of_coils==0) ){
		printf("\nWarning : 'upload_csm' : plan not initialized or number of coils is 0. Returning.\n");
		return false;
	}

	// Allocate device memory
	if( plan->CSMdevPtr )
		cudaFree( plan->CSMdevPtr );
	cudaMalloc( (void**) &plan->CSMdevPtr, plan->number_of_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex) );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'upload_csm': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Upload CSM
	cudaMemcpy( plan->CSMdevPtr, csmHostPtr, plan->number_of_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex), cudaMemcpyHostToDevice );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'upload_csm': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Normalize CSM.
	normalize_csm( plan->number_of_coils, prod(plan->matrix_size), plan->CSMdevPtr );

	// Calculate intensity correction image (preconditioning)
	if( plan->intensity_correction_magnitudes_image_DevPtr == 0x0 ){

		cudaMalloc( (void**) &plan->intensity_correction_magnitudes_image_DevPtr, prod(plan->matrix_size)*sizeof(float) );

		err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected in 'upload_csm': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}
	}

	calculate_intensity_correction_magnitudes_image( prod(plan->matrix_size), plan->number_of_coils, plan->CSMdevPtr, plan->intensity_correction_magnitudes_image_DevPtr );

/*
	//DEBUG
	float* _tmp = (float*) calloc( 1, plan->number_of_coils*prod(plan->matrix_size)*sizeof(float) );
	image_modulus( plan->CSMdevPtr, (float*) plan->CSMdevPtr, plan->number_of_coils*prod(plan->matrix_size), false, 1.0f );
	cudaMemcpy( _tmp, plan->CSMdevPtr, plan->number_of_coils*prod(plan->matrix_size)*sizeof(float), cudaMemcpyDeviceToHost );
	FILE *_fout = fopen("csm.raw", "wb");
	fwrite( _tmp, plan->number_of_coils*prod(plan->matrix_size), sizeof(float), _fout );
	fclose(_fout);
*/
	return true;
}

// Upload regularization/training data from the host to the device
template< class UINTd, class FLOATd, char TYPE > __host__ bool 
upload_regularization_image( mr_recon::NFFT_iteration_plan<UINTd,FLOATd, TYPE> *plan, float *image_HostPtr )
{
  if( !plan->initialized ){
    printf("\nInitilize the plan before calling 'upload_regularization_image'. Ignoring upload.\n");
    return false;
  }
  
  // Cleanup previous transfers
  if( plan->regularizationDevPtr ){
   cudaFree( plan->regularizationDevPtr );
   plan->regularizationDevPtr = 0x0;
  }
  
  cudaMalloc( (void**) &plan->regularizationDevPtr, prod(plan->matrix_size)*sizeof(float) );
  
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    printf("\nCuda error detected in 'upload_regularization_image': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
    exit(1);
  }
  
  // Copy host data to device
  cudaMemcpy( plan->regularizationDevPtr, image_HostPtr, prod(plan->matrix_size)*sizeof(float), cudaMemcpyHostToDevice );
  
  err = cudaGetLastError();
  if( err != cudaSuccess ){
    printf("\nCuda error detected in 'upload_regularization_image': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
    exit(1);
  }
  
  // Multiply instead of divide with regularization image
  image_reciprocal( prod(plan->matrix_size), plan->regularizationDevPtr );
  
  return true;
}


__global__ void 
rotate_csm_kernel( uint2 matrix_size, unsigned int num_coils, cuFloatComplex *in, cuFloatComplex *out )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int coil = blockIdx.y;

	if( idx < prod(matrix_size) ){
		
		uint2 out_co = idx_to_co(idx, matrix_size);
		const uint2 in_co = make_uint2( out_co.y, out_co.x );
		out[coil*prod(matrix_size)+idx] = in[coil*prod(matrix_size)+co_to_idx(in_co,matrix_size)];
	}
}

// TODO: templetize...
void
rotate_csm( mr_recon::NFFT_iteration_plan<uint2, float2, 2> *plan )
{
	if( !plan->CSMdevPtr )
	{
		printf("\nError in 'rotate_csm': No CSM computed or uploaded. Quitting.\n" ); fflush(stdout);	
		exit(1);
	}

	if( plan->matrix_size.x != plan->matrix_size.y )
	{
		printf("\nError in 'rotate_csm': only square matrices allowed. Quitting.\n" ); fflush(stdout);	
		exit(1);
	}

	// Make temporary copy of CSM
	cuFloatComplex *new_csm;
	cudaMalloc( (void**) &new_csm, plan->number_of_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex) );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'rotate_csm': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
	
	// Find dimensions of grid/blocks.
	dim3 dimBlock( 256 );
	dim3 dimGrid( (unsigned int) ceil((double) (prod(plan->matrix_size))/dimBlock.x), plan->number_of_coils );

	// Invoke kernel
	rotate_csm_kernel<<< dimGrid, dimBlock >>>( plan->matrix_size, plan->number_of_coils, plan->CSMdevPtr, new_csm );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'rotate_csm': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	cudaFree(plan->CSMdevPtr);
	plan->CSMdevPtr = new_csm;

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'rotate_csm': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	if( plan->intensity_correction_magnitudes_image_DevPtr == 0x0 ){

		cudaMalloc( (void**) &plan->intensity_correction_magnitudes_image_DevPtr, prod(plan->matrix_size)*sizeof(float) );

		err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected in 'upload_csm': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}
	}

	calculate_intensity_correction_magnitudes_image( prod(plan->matrix_size), plan->number_of_coils, plan->CSMdevPtr, plan->intensity_correction_magnitudes_image_DevPtr );

/*
	//DEBUG
	float* _tmp = (float*) calloc( 1, plan->number_of_coils*prod(plan->matrix_size)*sizeof(float) );
	image_modulus( plan->CSMdevPtr, (float*) plan->CSMdevPtr, plan->number_of_coils*prod(plan->matrix_size), false, 1.0f );
	cudaMemcpy( _tmp, plan->CSMdevPtr, plan->number_of_coils*prod(plan->matrix_size)*sizeof(float), cudaMemcpyDeviceToHost );
	FILE *_fout = fopen("csm.raw", "wb");
	fwrite( _tmp, plan->number_of_coils*prod(plan->matrix_size), sizeof(float), _fout );
	fclose(_fout);
*/
}

template< class UINTd, class FLOATd, char TYPE > bool
extract_csm_and_regularization( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE> *plan, float sigma_csm, unsigned int num_samples_per_coil, const cuFloatComplex *samples_HostPtr, float noise_decorrelation_shutter_radius, unsigned int global_projection_offset, unsigned int frames_per_rotation )
{
	// Some utility variables

	UINTd matrix_size = plan->matrix_size;
	UINTd matrix_size_os = plan->matrix_size_os;
	unsigned int num_griddings = num_samples_per_coil / plan->number_of_samples;

	/*
		Some initial checks for input/plan consistency
	*/

	if( sum(plan->fixed_dims) ){
		printf("\nWarning : 'extract_csm_and_regularization' : There can be no fixed dims. Returning.\n");
		return false;
	}

	if( !plan->initialized || (plan->number_of_coils==0) ){
		printf("\nWarning : 'extract_csm_and_regularization' : plan not initialized or number of coils is 0. Returning.\n");
		return false;
	}

	if( num_griddings == 0 ){
		printf("\nWarning : 'extract_csm_and_regularization' : There is not enough samples for the desired computation. Returning.\n");
		return false;	
	}

	if( num_griddings*plan->number_of_samples != num_samples_per_coil ){
		printf("\nWarning : 'extract_csm_and_regularization' : The number of samples is not a multiplum of the plan's required number of samples. The remainder will be ignored.\n");
	}

	// Be optimistic
	bool success = true;


	// These are the "return" data. Allocate.

	if( plan->CSMdevPtr )
		cudaFree( plan->CSMdevPtr );
	cudaMalloc( (void**) &plan->CSMdevPtr, plan->number_of_coils*prod(matrix_size)*sizeof(cuFloatComplex) );

	if( plan->regularizationDevPtr )
		cudaFree( plan->regularizationDevPtr );
	cudaMalloc( (void**) &plan->regularizationDevPtr, prod(matrix_size)*sizeof(float) );

	// Temporary storage
	cuFloatComplex *tmpRegularizationDevPtr;
	cudaMalloc( (void**) &tmpRegularizationDevPtr, plan->number_of_coils*prod(matrix_size)*sizeof(cuFloatComplex) );

	/*
		Compute coil sensitivity maps (and prepare for regularization computation)
	*/

	{
		// Temporary storage
		cuFloatComplex *accCSMdevPtr_os_DevPtr, *accRegularizationDevPtr_os_DevPtr, *samples_DevPtr;

		cudaMalloc((void**)&accCSMdevPtr_os_DevPtr, plan->number_of_coils*prod(matrix_size_os)*sizeof(cuFloatComplex));
		cudaMalloc((void**)&accRegularizationDevPtr_os_DevPtr, plan->number_of_coils*prod(matrix_size_os)*sizeof(cuFloatComplex));
		cudaMalloc((void**)&samples_DevPtr, plan->number_of_coils*plan->number_of_samples*sizeof(cuFloatComplex));

		// Error check
		cudaError_t err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected: %s\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}

		for( unsigned int iteration=0; iteration<num_griddings; iteration++ ){

			// Projection offset for current reconstruction
			unsigned int projection_offset = iteration*plan->projections_per_frame+global_projection_offset;

			// Run preprocessing for current reconstruction
			if( success )
				if( TYPE == 1 ) // golden angle
					success = preprocess_radial_NFFT( plan, plan->projections_per_frame, plan->samples_per_projection, plan->projections_per_frame, projection_offset );
				else
					success = preprocess_radial_NFFT( plan, plan->projections_per_frame, plan->samples_per_projection, plan->projections_per_frame, iteration%frames_per_rotation, frames_per_rotation );

			// Compute density compensation weights 
			if( success )
				success = compute_dcw_radial( plan );

			// Upload samples for reconstruction
			for( unsigned int c=0; c<plan->number_of_coils; c++ )
				cudaMemcpy( &samples_DevPtr[c*plan->number_of_samples], &samples_HostPtr[c*num_samples_per_coil+iteration*plan->number_of_samples], plan->number_of_samples*sizeof(cuFloatComplex), cudaMemcpyHostToDevice );

			//
			// !!! NOTE !!!
			// We assume that no noise decorrelation has been performed on the CPU
			// So we do it now...
			//
			noise_decorrelate_radial( plan, noise_decorrelation_shutter_radius, samples_DevPtr );

			// Density compensation
			density_compensate( plan, samples_DevPtr, samples_DevPtr );

			// Convolution
			for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 
				if( success )
					success = NFFT_iteration_convolve_to_image( plan, &samples_DevPtr[i*plan->domain_size_coils*plan->number_of_samples], &accCSMdevPtr_os_DevPtr[i*plan->domain_size_coils*prod(matrix_size_os)], (iteration>0) );
			}
		}

		// Take the average (?required for correctly scaled regularization?, not the CSM due to subsequent normalization)
//		cublasCscal( plan->number_of_coils*prod(matrix_size_os), make_cuFloatComplex(1.0f/num_griddings, 1.0f/num_griddings), accCSMdevPtr_os_DevPtr, 1 );

		// Copy to regularization temporary before smoothing
		cudaMemcpy( accRegularizationDevPtr_os_DevPtr, accCSMdevPtr_os_DevPtr, plan->number_of_coils*prod(plan->matrix_size_os)*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

		// Smoothen CSM (Gaussian convolution in image space)
		image_multiply_gaussian( matrix_size_os, sigma_csm*plan->alpha, accCSMdevPtr_os_DevPtr, plan->number_of_coils );

		// FFT accCSMdevPtr
		if( success )
			success = K2I_ALL( accCSMdevPtr_os_DevPtr, matrix_size_os, plan->number_of_coils, true, true );

		// FFT accRegularizationDevPtr_os_DevPtr
		if( success )
			success = K2I_ALL( accRegularizationDevPtr_os_DevPtr, matrix_size_os, plan->number_of_coils, true, true );

		// Crop
		if( success )
			crop_image( matrix_size, matrix_size_os, plan->CSMdevPtr, accCSMdevPtr_os_DevPtr, plan->number_of_coils );
		if( success )
			crop_image( matrix_size, matrix_size_os, tmpRegularizationDevPtr, accRegularizationDevPtr_os_DevPtr, plan->number_of_coils );

		// Deapodize
		for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 
			if( success )
				deapodize( plan, false, &plan->CSMdevPtr[i*plan->domain_size_coils*prod(matrix_size)], true );
			if( success )
				deapodize( plan, false, &tmpRegularizationDevPtr[i*plan->domain_size_coils*prod(matrix_size)] );
		}

		// Normalize and we have our CSM.
		normalize_csm( plan->number_of_coils, prod(matrix_size), plan->CSMdevPtr );

		cudaFree(samples_DevPtr); samples_DevPtr = 0x0;
		cudaFree(accRegularizationDevPtr_os_DevPtr); accRegularizationDevPtr_os_DevPtr = 0x0;
		cudaFree(accCSMdevPtr_os_DevPtr); accCSMdevPtr_os_DevPtr = 0x0;
	}

	/*
		Compute regularization
	*/

	// Initial complex estimate (Sense combination of the coils to a single image)
	extract_training_data_xt( plan->number_of_coils, prod(matrix_size), tmpRegularizationDevPtr, plan->CSMdevPtr, tmpRegularizationDevPtr );

	// Image modulus
	image_modulus( tmpRegularizationDevPtr, plan->regularizationDevPtr, prod(matrix_size), false, 1.0f );

	// Normalize to the square root of the number of image elements
	const float scale = sqrtf((float)prod(matrix_size))/cublasSasum( prod(matrix_size), plan->regularizationDevPtr, 1 );
	cublasSscal( prod(matrix_size), scale, plan->regularizationDevPtr, 1 );

	// Square
	square_image( prod(matrix_size), plan->regularizationDevPtr );

	// Multiply instead of divide with regularization image
	image_reciprocal( prod(matrix_size), plan->regularizationDevPtr );

/*
	//DEBUG
	float* _tmp = (float*) calloc( 1, prod(matrix_size)*sizeof(float) );
	cudaMemcpy( _tmp, plan->regularizationDevPtr, prod(matrix_size)*sizeof(float), cudaMemcpyDeviceToHost );
	FILE *_fout = fopen("theta.raw", "wb");
	fwrite( _tmp, prod(matrix_size), sizeof(float), _fout );
	fclose(_fout);
*/

	// Calculate intensity correction image (preconditioning)
	if( plan->intensity_correction_magnitudes_image_DevPtr == 0x0 )
		cudaMalloc( (void**) &plan->intensity_correction_magnitudes_image_DevPtr, prod(matrix_size)*sizeof(float) );	
	calculate_intensity_correction_magnitudes_image( prod(matrix_size), plan->number_of_coils, plan->CSMdevPtr, plan->intensity_correction_magnitudes_image_DevPtr );

	// Cleanup
	cudaFree(tmpRegularizationDevPtr); tmpRegularizationDevPtr = 0x0;

/*
	//DEBUG
	float* __tmp = (float*) calloc( 1, plan->number_of_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex) );
	cudaMemcpy( __tmp, plan->CSMdevPtr, plan->number_of_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("csm.raw", "wb");
	fwrite( __tmp, plan->number_of_coils*prod(plan->matrix_size), sizeof(cuFloatComplex), fout );
	fclose(fout);
*/

	return true;
}


__global__ void
rms_combine_kernel( unsigned int num_coils, unsigned int num_elements, cuFloatComplex *in, cuFloatComplex* out)
{

	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if( idx < num_elements ){
		float sum_sq = 0.0f;

		for( unsigned int c=0; c<num_coils; c++){
			cuFloatComplex val = in[c*num_elements+idx];
			sum_sq += val.x*val.x+val.y*val.y;
		}


		out[idx] = make_cuFloatComplex(sqrt(sum_sq),0.0f);
	}
}

void
rms_combine( unsigned int num_coils, unsigned int num_elements, cuFloatComplex *in, cuFloatComplex *out )
{
	dim3 blockDim(512,1,1);
	dim3 gridDim((unsigned int) ceil((double)(num_elements/blockDim.x)), 1, 1 );

	assert(num_elements>=blockDim.x);

	rms_combine_kernel<<< gridDim, blockDim >>>( num_coils, num_elements, in, out );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'normalize_csm_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

template< class UINTd, class FLOATd, char TYPE > bool
extract_and_combine_current_frame( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE> *plan, cuFloatComplex *kspace_coil_images_buffer_os,  cuFloatComplex *imageDevPtr)
{
	bool success = true;

	UINTd matrix_size = plan->matrix_size;
	UINTd matrix_size_os = plan->matrix_size_os;
	unsigned int number_of_coils = plan->number_of_coils;

      	cuFloatComplex* tmp_grid_os;
      	cuFloatComplex* tmp_grid;
      	cudaMalloc( (void**) &tmp_grid_os, number_of_coils*prod(matrix_size_os)*sizeof(cuFloatComplex) );
      	cudaMalloc( (void**) &tmp_grid, number_of_coils*prod(matrix_size)*sizeof(cuFloatComplex) );

	cudaError_t err = cudaGetLastError();
      	if( err != cudaSuccess ){
	    printf("extract_and_combine_current_frame: Cuda error detected in memory allocation (temprary gridding memory): %s\n", cudaGetErrorString(err));
	    fflush(stdout);
	    exit(1);
      	}


        cudaMemcpy(tmp_grid_os, kspace_coil_images_buffer_os, prod(matrix_size_os)*number_of_coils*sizeof(cuFloatComplex),cudaMemcpyDeviceToDevice); 

        success = K2I_ALL( tmp_grid_os, matrix_size_os, number_of_coils, true, true );

        if( success )
	    crop_image( matrix_size, matrix_size_os, tmp_grid, tmp_grid_os, number_of_coils );

	// Deapodize
	for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 
		if( success )
		    deapodize( plan, false, &tmp_grid[i*plan->domain_size_coils*prod(matrix_size)], true );
	}


	rms_combine( number_of_coils, prod(matrix_size), tmp_grid, imageDevPtr );

      
	cudaFree(tmp_grid_os);
      	cudaFree(tmp_grid);


	return success;

}

template< class UINTd, class FLOATd, char TYPE > bool
update_csm_and_regularization( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE> *plan, const cuFloatComplex *_samples_DevPtr, float shutter_csm, cuFloatComplex *kspace_coil_images_buffer_os, unsigned int projection_offset, bool do_preprocessing, unsigned int frames_per_rotation, bool reallocate )
{
	// Some utility variables
	UINTd matrix_size = plan->matrix_size;
	UINTd matrix_size_os = plan->matrix_size_os;

	// Some initial checks for input/plan consistency

	if( sum(plan->fixed_dims) ){
		printf("\nWarning : 'update_csm_and_regularization' : There can be no fixed dims. Returning.\n");
		return false;
	}

	if( !plan->initialized || (plan->number_of_coils==0) ){
		printf("\nWarning : 'update_csm_and_regularization' : plan not initialized or number of coils is 0. Returning.\n");
		return false;
	}

	// These are the "return" data. Allocate.

	if( plan->CSMdevPtr && reallocate ){
		cudaFree( plan->CSMdevPtr );
		plan->CSMdevPtr = 0x0;
	}

	if( plan->CSMdevPtr == 0x0 )
		cudaMalloc( (void**) &plan->CSMdevPtr, plan->number_of_coils*prod(matrix_size)*sizeof(cuFloatComplex) );

	if( plan->regularizationDevPtr && reallocate ){
		cudaFree( plan->regularizationDevPtr );
		plan->regularizationDevPtr = 0x0;
	}

	if( plan->regularizationDevPtr == 0x0 )
		cudaMalloc( (void**) &plan->regularizationDevPtr, prod(matrix_size)*sizeof(float) );

	// Static storage
	static cuFloatComplex *kspace_acc_coil_images_os[] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
	static unsigned int kspace_acc_coil_images_os_size[] = {0,0,0,0,0,0,0,0};
	
	int dev_number;
	cudaGetDevice(&dev_number);

	if( kspace_acc_coil_images_os[dev_number] == 0x0 ){
		cudaMalloc( (void**) &kspace_acc_coil_images_os[dev_number], prod(matrix_size_os)*plan->number_of_coils*sizeof(cuFloatComplex) );
		clear_image( prod(matrix_size_os)*plan->number_of_coils, make_cuFloatComplex(0.0f, 0.0f), kspace_acc_coil_images_os[dev_number] );
		kspace_acc_coil_images_os_size[dev_number] = prod(matrix_size_os)*plan->number_of_coils;
	}

	// If 'matrix_size_os' or the number of coils change we need to update 'kspace_acc_coil_images_os'
	if( kspace_acc_coil_images_os_size[dev_number] != prod(matrix_size_os)*plan->number_of_coils ){
		if( kspace_acc_coil_images_os[dev_number] )
			cudaFree(kspace_acc_coil_images_os[dev_number]);
		cudaMalloc( (void**) &kspace_acc_coil_images_os[dev_number], prod(matrix_size_os)*plan->number_of_coils*sizeof(cuFloatComplex) );
		clear_image( prod(matrix_size_os)*plan->number_of_coils, make_cuFloatComplex(0.0f, 0.0f), kspace_acc_coil_images_os[dev_number] );
		kspace_acc_coil_images_os_size[dev_number] = prod(matrix_size_os)*plan->number_of_coils;
	}

	// Temporary storage
	cuFloatComplex *tmpRegularizationDevPtr = 0x0;
	cudaMalloc( (void**) &tmpRegularizationDevPtr, plan->number_of_coils*prod(matrix_size)*sizeof(cuFloatComplex) );

	// Error check
	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected: 'update_csm_and_regularization' : %s\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Subtract the buffer we are about to replace from the accumulated images
	subtract_images( prod(matrix_size_os)*plan->number_of_coils, kspace_acc_coil_images_os[dev_number], kspace_acc_coil_images_os[dev_number], kspace_coil_images_buffer_os );

	// Temporary storage
	cuFloatComplex *tmp_acc_os_DevPtr = 0x0, *samples_DevPtr = 0x0;
	cudaMalloc((void**)&tmp_acc_os_DevPtr, plan->number_of_coils*prod(matrix_size_os)*sizeof(cuFloatComplex));
	cudaMalloc((void**)&samples_DevPtr, plan->number_of_coils*plan->number_of_samples*sizeof(cuFloatComplex));

	// Error check
	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected: 'update_csm_and_regularization' : %s\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Copy samples before density compensation
	cudaMemcpy( samples_DevPtr, _samples_DevPtr, plan->number_of_samples*plan->number_of_coils*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

	// Be optimistic
	bool success = true;

	/*
		Compute contribution from current frame
	*/

	// Run preprocessing for current reconstruction
	if( success && do_preprocessing )
		success = preprocess_radial_NFFT( plan, plan->projections_per_frame, plan->samples_per_projection, plan->projections_per_frame, projection_offset, frames_per_rotation );

	// Compute density compensation weights 
	if( success )
		success = compute_dcw_radial( plan );

	// Density compensation
	if( success )
		success = density_compensate( plan, samples_DevPtr, samples_DevPtr );

	// Convolutions
	for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 
		if( success )
			success = NFFT_iteration_convolve_to_image( plan, &samples_DevPtr[i*plan->domain_size_coils*plan->number_of_samples], &kspace_coil_images_buffer_os[i*plan->domain_size_coils*prod(matrix_size_os)] );
	}

	// Add buffers we just computed to the accumulated image
	add_images( prod(matrix_size_os)*plan->number_of_coils, kspace_acc_coil_images_os[dev_number], kspace_acc_coil_images_os[dev_number], kspace_coil_images_buffer_os );

	/*
		Compute CSM
	*/

	// Make temporary image before smoothing
	cudaMemcpy( tmp_acc_os_DevPtr, kspace_acc_coil_images_os[dev_number], plan->number_of_coils*prod(plan->matrix_size_os)*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

	// Smoothen CSM (Gaussian convolution in image space)
	image_multiply_gaussian( matrix_size_os, shutter_csm*plan->alpha, tmp_acc_os_DevPtr, plan->number_of_coils );

	// Error check
	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected: 'update_csm_and_regularization' : %s\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// FFT CSM
	if( success )
		success = K2I_ALL( tmp_acc_os_DevPtr, matrix_size_os, plan->number_of_coils, true, true );

	// Crop
	if( success )
		crop_image( matrix_size, matrix_size_os, plan->CSMdevPtr, tmp_acc_os_DevPtr, plan->number_of_coils );

	// Copy acc buffer to regularization temporary before FFT. 
	// We really ought to make an FFT wrapper for the "out of place" case as well...
	cudaMemcpy( tmp_acc_os_DevPtr, kspace_acc_coil_images_os[dev_number], plan->number_of_coils*prod(plan->matrix_size_os)*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

	// FFT
	if( success )
		success = K2I_ALL( tmp_acc_os_DevPtr, matrix_size_os, plan->number_of_coils, true, true );

	if( success )
		crop_image( matrix_size, matrix_size_os, tmpRegularizationDevPtr, tmp_acc_os_DevPtr, plan->number_of_coils );

	// Deapodize
	for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 
		if( success )
			deapodize( plan, false, &plan->CSMdevPtr[i*plan->domain_size_coils*prod(matrix_size)], true );
		if( success )
			deapodize( plan, false, &tmpRegularizationDevPtr[i*plan->domain_size_coils*prod(matrix_size)] );
	}

	// Normalize and we have our CSM.
	normalize_csm( plan->number_of_coils, prod(matrix_size), plan->CSMdevPtr );

	cudaFree(samples_DevPtr); samples_DevPtr = 0x0;
	cudaFree(tmp_acc_os_DevPtr); tmp_acc_os_DevPtr = 0x0;

	/*
		Compute regularization
	*/

	// Initial complex estimate (SENSE combination of the coils to a single image)
	extract_training_data_xt( plan->number_of_coils, prod(matrix_size), tmpRegularizationDevPtr, plan->CSMdevPtr, tmpRegularizationDevPtr );

	// Image modulus
	image_modulus( tmpRegularizationDevPtr, plan->regularizationDevPtr, prod(matrix_size), false, 1.0f );

	// Normalize to the square root of the number of image elements
	const float scale = sqrtf((float)prod(matrix_size))/cublasSasum( prod(matrix_size), plan->regularizationDevPtr, 1 );
	cublasSscal( prod(matrix_size), scale, plan->regularizationDevPtr, 1 );

	// Squared
	square_image( prod(matrix_size), plan->regularizationDevPtr );

	//clear_image( prod(matrix_size), 1.0f, plan->regularizationDevPtr );

	/*
	//DEBUG
	float* _tmp = (float*) calloc( 1, prod(matrix_size)*sizeof(float) );
	cudaMemcpy( _tmp, plan->regularizationDevPtr, prod(matrix_size)*sizeof(float), cudaMemcpyDeviceToHost );
	FILE *_fout = fopen("theta.raw", "wb");
	fwrite( _tmp, prod(matrix_size), sizeof(float), _fout );
	fclose(_fout);
	*/

	// Multiply instead of divide with regularization image
	image_reciprocal( prod(matrix_size), plan->regularizationDevPtr );

	// Calculate intensity correction image (preconditioning)
	if( plan->intensity_correction_magnitudes_image_DevPtr == 0x0 )
		cudaMalloc( (void**) &plan->intensity_correction_magnitudes_image_DevPtr, prod(matrix_size)*sizeof(float) );	
	calculate_intensity_correction_magnitudes_image( prod(matrix_size), plan->number_of_coils, plan->CSMdevPtr, plan->intensity_correction_magnitudes_image_DevPtr );

	// Cleanup
	cudaFree(tmpRegularizationDevPtr); tmpRegularizationDevPtr = 0x0;

	/*
	//DEBUG
	float* __tmp = (float*) calloc( 1, plan->number_of_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex) );
	cudaMemcpy( __tmp, plan->CSMdevPtr, plan->number_of_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("csm.raw", "wb");
	fwrite( __tmp, plan->number_of_coils*prod(plan->matrix_size), sizeof(cuFloatComplex), fout );
	fclose(fout);
	*/

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected: 'update_csm_and_regularization' : %s\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	return success;
}

template< class UINTd, class FLOATd > bool
update_csm_and_regularization( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, 0> *plan, cuFloatComplex *_samples_DevPtr, FLOATd *traj_DevPtr, float *dcw_DevPtr, float shutter_csm, cuFloatComplex *kspace_coil_images_buffer_os, cuFloatComplex *kspace_acc_coil_images_os, unsigned int kspace_acc_coil_images_os_size, bool do_preprocessing, bool reallocate )
{
	// Some utility variables
	UINTd matrix_size = plan->matrix_size;
	UINTd matrix_size_os = plan->matrix_size_os;

	// Some initial checks for input/plan consistency

	if( sum(plan->fixed_dims) ){
		printf("\nWarning : 'update_csm_and_regularization' : There can be no fixed dims. Returning.\n");
		return false;
	}

	if( !plan->initialized || (plan->number_of_coils==0) ){
		printf("\nWarning : 'update_csm_and_regularization' : plan not initialized or number of coils is 0. Returning.\n");
		return false;
	}

	// These are the "return" data. Allocate.

	if( plan->CSMdevPtr && reallocate ){
		cudaFree( plan->CSMdevPtr );
		plan->CSMdevPtr = 0x0;
	}

	if( plan->CSMdevPtr == 0x0 )
		cudaMalloc( (void**) &plan->CSMdevPtr, plan->number_of_coils*prod(matrix_size)*sizeof(cuFloatComplex) );

	if( plan->regularizationDevPtr && reallocate ){
		cudaFree( plan->regularizationDevPtr );
		plan->regularizationDevPtr = 0x0;
	}

	if( plan->regularizationDevPtr == 0x0 )
		cudaMalloc( (void**) &plan->regularizationDevPtr, prod(matrix_size)*sizeof(float) );

	// Static storage
	/*
	static cuFloatComplex *kspace_acc_coil_images_os[] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
	static unsigned int kspace_acc_coil_images_os_size[] = {0,0,0,0,0,0,0,0};

	int dev_number;
	cudaGetDevice(&dev_number);

	if( kspace_acc_coil_images_os[dev_number] == 0x0 ){
		cudaMalloc( (void**) &kspace_acc_coil_images_os[dev_number], prod(matrix_size_os)*plan->number_of_coils*sizeof(cuFloatComplex) );
		clear_image( prod(matrix_size_os)*plan->number_of_coils, make_cuFloatComplex(0.0f, 0.0f), kspace_acc_coil_images_os[dev_number] );
		kspace_acc_coil_images_os_size[dev_number] = prod(matrix_size_os)*plan->number_of_coils;
	}

	// If 'matrix_size_os' or the number of coils change we need to update 'kspace_acc_coil_images_os'
	if( kspace_acc_coil_images_os_size[dev_number] != prod(matrix_size_os)*plan->number_of_coils ){
		if( kspace_acc_coil_images_os[dev_number] )
			cudaFree(kspace_acc_coil_images_os[dev_number]);
		cudaMalloc( (void**) &kspace_acc_coil_images_os[dev_number], prod(matrix_size_os)*plan->number_of_coils*sizeof(cuFloatComplex) );
		clear_image( prod(matrix_size_os)*plan->number_of_coils, make_cuFloatComplex(0.0f, 0.0f), kspace_acc_coil_images_os[dev_number] );
		kspace_acc_coil_images_os_size[dev_number] = prod(matrix_size_os)*plan->number_of_coils;
	}

	*/


	// Temporary storage
	cuFloatComplex *tmpRegularizationDevPtr = 0x0;
	cudaMalloc( (void**) &tmpRegularizationDevPtr, plan->number_of_coils*prod(matrix_size)*sizeof(cuFloatComplex) );

	// Error check
	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected: 'update_csm_and_regularization' : %s\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Subtract the buffer we are about to replace from the accumulated images
	//subtract_images( prod(matrix_size_os)*plan->number_of_coils, kspace_acc_coil_images_os[dev_number], kspace_acc_coil_images_os[dev_number], kspace_coil_images_buffer_os );
	subtract_images( prod(matrix_size_os)*plan->number_of_coils, kspace_acc_coil_images_os, kspace_acc_coil_images_os, kspace_coil_images_buffer_os );

	// Temporary storage
	cuFloatComplex *tmp_acc_os_DevPtr = 0x0, *samples_DevPtr = 0x0;
	cudaMalloc((void**)&tmp_acc_os_DevPtr, plan->number_of_coils*prod(matrix_size_os)*sizeof(cuFloatComplex));
	cudaMalloc((void**)&samples_DevPtr, plan->number_of_coils*plan->number_of_samples*sizeof(cuFloatComplex));

	// Error check
	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected: 'update_csm_and_regularization' : %s\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Copy samples before density compensation
	cudaMemcpy( samples_DevPtr, _samples_DevPtr, plan->number_of_samples*plan->number_of_coils*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

	// Be optimistic
	bool success = true;

	/*
		Compute contribution from current frame
	*/

	// Run preprocessing for current reconstruction
	if( success && do_preprocessing )
		success = preprocess_generic_NFFT( plan, plan->number_of_samples, traj_DevPtr );

	// Compute density compensation weights 
	if( success )
		set_dcw( plan, dcw_DevPtr );

	// Density compensation
	if( success )
		success = density_compensate( plan, samples_DevPtr, samples_DevPtr );

	// Convolutions
	for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 
		if( success )
			success = NFFT_iteration_convolve_to_image( plan, &samples_DevPtr[i*plan->domain_size_coils*plan->number_of_samples], &kspace_coil_images_buffer_os[i*plan->domain_size_coils*prod(matrix_size_os)] );
	}

	// Add buffers we just computed to the accumulated image
	add_images( prod(matrix_size_os)*plan->number_of_coils, kspace_acc_coil_images_os, kspace_acc_coil_images_os, kspace_coil_images_buffer_os );

	/*
		Compute CSM
	*/

	// Make temporary image before smoothing
	cudaMemcpy( tmp_acc_os_DevPtr, kspace_acc_coil_images_os, plan->number_of_coils*prod(plan->matrix_size_os)*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

	// Smoothen CSM (Gaussian convolution in image space)
	image_multiply_gaussian( matrix_size_os, shutter_csm*plan->alpha, tmp_acc_os_DevPtr, plan->number_of_coils );

	// Error check
	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected: 'update_csm_and_regularization' : %s\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// FFT CSM
	if( success )
		success = K2I_ALL( tmp_acc_os_DevPtr, matrix_size_os, plan->number_of_coils, true, true );

	// Crop
	if( success )
		crop_image( matrix_size, matrix_size_os, plan->CSMdevPtr, tmp_acc_os_DevPtr, plan->number_of_coils );

	// Copy acc buffer to regularization temporary before FFT. 
	// We really ought to make an FFT wrapper for the "out of place" case as well...
	cudaMemcpy( tmp_acc_os_DevPtr, kspace_acc_coil_images_os, plan->number_of_coils*prod(plan->matrix_size_os)*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

	// FFT
	if( success )
		success = K2I_ALL( tmp_acc_os_DevPtr, matrix_size_os, plan->number_of_coils, true, true );

	if( success )
		crop_image( matrix_size, matrix_size_os, tmpRegularizationDevPtr, tmp_acc_os_DevPtr, plan->number_of_coils );

	// Deapodize
	for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 
		if( success )
			deapodize( plan, false, &plan->CSMdevPtr[i*plan->domain_size_coils*prod(matrix_size)], true );
		if( success )
			deapodize( plan, false, &tmpRegularizationDevPtr[i*plan->domain_size_coils*prod(matrix_size)] );
	}

	// Normalize and we have our CSM.
	normalize_csm( plan->number_of_coils, prod(matrix_size), plan->CSMdevPtr );

	cudaFree(samples_DevPtr); samples_DevPtr = 0x0;
	cudaFree(tmp_acc_os_DevPtr); tmp_acc_os_DevPtr = 0x0;

	/*
		Compute regularization
	*/

	// Initial complex estimate (SENSE combination of the coils to a single image)
	extract_training_data_xt( plan->number_of_coils, prod(matrix_size), tmpRegularizationDevPtr, plan->CSMdevPtr, tmpRegularizationDevPtr );

	// Image modulus
	image_modulus( tmpRegularizationDevPtr, plan->regularizationDevPtr, prod(matrix_size), false, 1.0f );

	// Normalize to the square root of the number of image elements
	const float scale = sqrtf((float)prod(matrix_size))/cublasSasum( prod(matrix_size), plan->regularizationDevPtr, 1 );
	cublasSscal( prod(matrix_size), scale, plan->regularizationDevPtr, 1 );

	// Squared
	square_image( prod(matrix_size), plan->regularizationDevPtr );

	//clear_image( prod(matrix_size), 1.0f, plan->regularizationDevPtr );

	/*
	//DEBUG
	float* _tmp = (float*) calloc( 1, prod(matrix_size)*sizeof(float) );
	cudaMemcpy( _tmp, plan->regularizationDevPtr, prod(matrix_size)*sizeof(float), cudaMemcpyDeviceToHost );
	FILE *_fout = fopen("theta.raw", "wb");
	fwrite( _tmp, prod(matrix_size), sizeof(float), _fout );
	fclose(_fout);
	*/

	// Multiply instead of divide with regularization image
	image_reciprocal( prod(matrix_size), plan->regularizationDevPtr );

	// Calculate intensity correction image (preconditioning)
	if( plan->intensity_correction_magnitudes_image_DevPtr == 0x0 )
		cudaMalloc( (void**) &plan->intensity_correction_magnitudes_image_DevPtr, prod(matrix_size)*sizeof(float) );	
	calculate_intensity_correction_magnitudes_image( prod(matrix_size), plan->number_of_coils, plan->CSMdevPtr, plan->intensity_correction_magnitudes_image_DevPtr );

	// Cleanup
	cudaFree(tmpRegularizationDevPtr); tmpRegularizationDevPtr = 0x0;

	/*
	//DEBUG
	float* __tmp = (float*) calloc( 1, plan->number_of_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex) );
	cudaMemcpy( __tmp, plan->CSMdevPtr, plan->number_of_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("csm.raw", "wb");
	fwrite( __tmp, plan->number_of_coils*prod(plan->matrix_size), sizeof(cuFloatComplex), fout );
	fclose(fout);
	*/

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected: 'update_csm_and_regularization' : %s\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	return success;
}


template< class UINTd, class FLOATd, char TYPE > bool
extract_csm_and_training_data( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE> *plan, float sigma_csm, unsigned int csm_frames, const cuFloatComplex *_samplesDevPtr, float _shutter_theta )
{
	// Some utility variables

	UINTd matrix_size = plan->matrix_size;
	UINTd matrix_size_os = plan->matrix_size_os;

	/*
		Some initial checks for input/plan consistency
	*/

	if( sum(plan->fixed_dims) != 1 || !get_last_dim(plan->fixed_dims) ){
		printf("\nWarning : 'extract_csm_and_training_data' : plan's last dimension must be fixed. There cannot be other fixed dimensions. Returning.");
		return false;
	}

	if( !plan->initialized || (plan->number_of_coils==0) ){
		printf("\nWarning : 'extract_csm_and_training_data' : plan not initialized or number of coils is 0. Returning.");
		return false;
	}

	if( !plan->weights_DevPtr ){
		printf("\nWarning : 'extract_csm_and_training_data' : density compensation weights not uploaded or computed. Returning.");
		return false;	
	}

	if( csm_frames > get_last_dim(matrix_size) ){
		printf("\nWarning : 'extract_csm_and_training_data' : csm_frames cannot exceed last dimension of matrix_size. Returning.");
		return false;	
	}

	// Be optimistic
	bool success = true;


	// These are the "return" data

	if( !plan->CSMdevPtr )
		cudaMalloc( (void**) &plan->CSMdevPtr, plan->number_of_coils*prod(matrix_size)*sizeof(cuFloatComplex) );
	if( !plan->regularizationDevPtr )
		cudaMalloc( (void**) &plan->regularizationDevPtr, prod(matrix_size)*sizeof(float) );

	float shutter_theta = _shutter_theta;
	if( shutter_theta < 0.0f )
		shutter_theta = compute_training_data_shutter_radius<TYPE>( plan->projections_per_frame );

//	printf("\nTraining data shutter radius is set to %f", shutter_theta );

	// Extract...

	cuFloatComplex *lowres, *hires;
	cudaMalloc((void**)&lowres, plan->number_of_coils*prod(matrix_size)*sizeof(cuFloatComplex));
	cudaMalloc((void**)&hires, plan->number_of_coils*prod(crop_last_dim(matrix_size))*sizeof(cuFloatComplex));

	/*
		Compute unaliased low-resolution images for each coil.
	*/

	{
		cuFloatComplex *lowres_os_DevPtr, *hires_os_DevPtr, *samples_DevPtr;

		cudaMalloc((void**)&lowres_os_DevPtr, plan->number_of_coils*prod(matrix_size_os)*sizeof(cuFloatComplex));
		cudaMalloc((void**)&hires_os_DevPtr, plan->number_of_coils*prod(crop_last_dim(matrix_size_os))*sizeof(cuFloatComplex));
		cudaMalloc((void**)&samples_DevPtr, plan->number_of_coils*plan->number_of_samples*sizeof(cuFloatComplex));

		// Error check
		cudaError_t err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected: %s\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}

		// Density compensation
		density_compensate( plan, samples_DevPtr, _samplesDevPtr );

		// Convolution
		for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 
			if( success )
				success = NFFT_iteration_convolve_to_image( plan, &samples_DevPtr[i*plan->domain_size_coils*plan->number_of_samples], &lowres_os_DevPtr[i*plan->domain_size_coils*prod(matrix_size_os)]);
		}

		// Copy k-space to hi-res time averaged (accumulated)
		for( unsigned int i=0; i<plan->number_of_coils; i++ )
			add_images( prod(crop_last_dim(matrix_size_os)), &hires_os_DevPtr[i*prod(crop_last_dim(matrix_size_os))], &lowres_os_DevPtr[i*get_last_dim(matrix_size_os)*prod(crop_last_dim(matrix_size_os))], csm_frames );

		// Circular shutter lowres (preserve only the fully sampled central part of kspace)
		cuda_border_fill_image<cuFloatComplex, uint2>( crop_last_dim(matrix_size_os), shutter_theta*plan->alpha, make_cuFloatComplex(0.0f, 0.0f), lowres_os_DevPtr, get_last_dim(matrix_size)*plan->number_of_coils );

		// Smoothen CSM (Gaussian convolution in image space)
		image_multiply_gaussian( crop_last_dim(matrix_size_os), sigma_csm*plan->alpha, hires_os_DevPtr, plan->number_of_coils );

		// FFT lowres (num_coils*num_images 2D transforms)
		if( success )
			success = K2I_ALL( lowres_os_DevPtr, crop_last_dim(matrix_size_os), plan->number_of_coils*get_last_dim(matrix_size_os), true, true );

		// Smoothing (in k-space) due to "steep" shutter cut-off (doesn't seem to really have an effect)
		//image_multiply_gaussian( crop_last_dim(matrix_size_os), ??? radius ???, lowres_os_DevPtr, plan->number_of_coils*get_last_dim(matrix_size_os) );

		// FFT hires (num_coils 2D transforms)
		if( success )
			success = K2I_ALL( hires_os_DevPtr, crop_last_dim(matrix_size_os), plan->number_of_coils, true, true );	

/*
	//DEBUG
	int num_elements = plan->number_of_coils*prod(matrix_size_os);
	float* _tmp = (float*) calloc( 1, num_elements*sizeof(cuFloatComplex) );
	cudaMemcpy( _tmp, lowres_os_DevPtr, num_elements*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("out.raw", "wb");
	fwrite( _tmp, num_elements, sizeof(cuFloatComplex), fout );
	fclose(fout);
*/

		// Crop
		if( success )
			crop_image( matrix_size, matrix_size_os, lowres, lowres_os_DevPtr, plan->number_of_coils );
		if( success )
			crop_image( crop_last_dim(matrix_size), crop_last_dim(matrix_size_os), hires, hires_os_DevPtr, plan->number_of_coils );

		// Deapodize
		for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 
			if( success )
				deapodize( plan, false, &lowres[i*plan->domain_size_coils*prod(matrix_size)] );
			if( success )
				deapodize( plan, false, &hires[i*plan->domain_size_coils*prod(crop_last_dim(matrix_size))], true );
		}

		// Normalize and we have our CSM.
		normalize_csm( plan->number_of_coils, prod(crop_last_dim(matrix_size)), hires );

		// Finally expand csm to the matrix_size (copy single csm to all frames for convenience).
		expand_csm( plan->number_of_coils, prod(crop_last_dim(matrix_size)), get_last_dim(matrix_size), hires, plan->CSMdevPtr );

		cudaFree(samples_DevPtr); samples_DevPtr = 0x0;
		cudaFree(lowres_os_DevPtr); lowres_os_DevPtr = 0x0;
		cudaFree(hires_os_DevPtr); hires_os_DevPtr = 0x0;
	}

	/*
		Compute training data (regularization)
	*/

	{		
		// Initial complex estimate in xt-space
		extract_training_data_xt( plan->number_of_coils, prod(matrix_size), lowres, plan->CSMdevPtr, lowres );

		// Transform training data to xf-space
		if( success )
			success = K2I( lowres, uintd_to_uint4_with_ones(matrix_size), 2, true );	

		// Image modulus
		image_modulus( lowres, plan->regularizationDevPtr, prod(matrix_size), false, 1.0f );
	
		// Normalize to the square root of the number of image elements
//		const float scale = fsqrtf((float)prod(matrix_size))/cublasSasum( prod(matrix_size), plan->regularizationDevPtr, 1 );
//		cublasSscal( prod(matrix_size), scale, plan->regularizationDevPtr, 1 );

		// Squared
		square_image( prod(matrix_size), plan->regularizationDevPtr );

/*
	//DEBUG
	float* _tmp = (float*) calloc( 1, prod(matrix_size)*sizeof(float) );
	cudaMemcpy( _tmp, regularizationDevPtr, prod(matrix_size)*sizeof(float), cudaMemcpyDeviceToHost );
	FILE *_fout = fopen("theta.raw", "wb");
	fwrite( _tmp, prod(matrix_size), sizeof(float), _fout );
	fclose(_fout);
*/

		// Multiply instead of divide with regularization image
		image_reciprocal( prod(matrix_size), plan->regularizationDevPtr );

		// Calculate intensity correction image
		if( plan->intensity_correction_magnitudes_image_DevPtr == 0x0 )
			cudaMalloc( (void**) &plan->intensity_correction_magnitudes_image_DevPtr, prod(matrix_size)*sizeof(float) );	
		calculate_intensity_correction_magnitudes_image( prod(matrix_size), plan->number_of_coils, plan->CSMdevPtr, plan->intensity_correction_magnitudes_image_DevPtr );

	}

	//DEBUG
/*
	float* __tmp = (float*) calloc( 1, num_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex) );
	for( unsigned int i=0; i<num_coils; i++ )
	cudaMemcpy( __tmp, CSMdevPtr, num_coils*prod(plan->matrix_size)*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("csm.raw", "wb");
	fwrite( __tmp, num_coils*prod(plan->matrix_size), sizeof(cuFloatComplex), fout );
	fclose(fout);
*/
	// Cleanup
	cudaFree(lowres); lowres = 0x0;
	cudaFree(hires); hires = 0x0;

	return true;
}

__global__ void
normalize_csm_kernel( unsigned int num_coils, unsigned int num_elements, cuFloatComplex *csm )
{

	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if( idx < num_elements ){
		
		cuFloatComplex sum_sq_cplx = make_cuFloatComplex(0.0f, 0.0f);;

		for( unsigned int c=0; c<num_coils; c++){
			cuFloatComplex val = csm[c*num_elements+idx];
			sum_sq_cplx = cuCaddf( sum_sq_cplx, cuCmulf( val, cuConjf(val)));
		}

		const float denominator = sqrtf(cuCrealf(sum_sq_cplx));
		float scale;
		
		if( denominator > 0.0001f )
			scale = 1.0f/denominator;
		else
			scale = 0.0f;

		for( unsigned int c=0; c<num_coils; c++){
			cuFloatComplex val = csm[c*num_elements+idx];
			csm[c*num_elements+idx] = make_cuFloatComplex( scale*cuCrealf(val), scale*cuCimagf(val) );
		}
	}
}

void
normalize_csm( unsigned int num_coils, unsigned int num_elements, cuFloatComplex *csm )
{
	dim3 blockDim(512,1,1);
	dim3 gridDim((unsigned int) ceil((double)(num_elements/blockDim.x)), 1, 1 );

	assert(num_elements>=blockDim.x);

	normalize_csm_kernel<<< gridDim, blockDim >>>( num_coils, num_elements, csm );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'normalize_csm_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

__global__ void 
expand_csm_kernel( unsigned int num_coils, unsigned int num_elements, unsigned int num_frames, const cuFloatComplex *in_csm, cuFloatComplex *out_csm )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if( idx < num_elements ){
		
		for( unsigned int c=0; c<num_coils; c++){	
			cuFloatComplex val = in_csm[c*num_elements+idx];
			for( unsigned int f=0; f<num_frames; f++)
				out_csm[c*num_frames*num_elements+f*num_elements+idx] = val;
		}
	}
}

void 
expand_csm( unsigned int num_coils, unsigned int num_elements, unsigned int num_frames, const cuFloatComplex *in_csm, cuFloatComplex *out_csm )
{
	dim3 blockDim(512,1,1);
	dim3 gridDim((unsigned int) ceil((double)(num_elements/blockDim.x)), 1, 1 );

	assert(num_elements>=blockDim.x);

	expand_csm_kernel<<< gridDim, blockDim >>>( num_coils, num_elements, num_frames, in_csm, out_csm );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'expand_csm_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


__global__ void 
extract_training_data_xt_kernel( unsigned int num_coils, unsigned int num_elements, const cuFloatComplex *lowres, const cuFloatComplex *csm, cuFloatComplex *ctraining )
{
	// This is basically a "fully sampled" non-Cartesian SENSE reconstruction.
	// I.e. the right hand side computation.
	// The norm ensures a "scaling factor" of 1.

	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if( idx < num_elements ){
		
		cuFloatComplex norm = make_cuFloatComplex(0.0f, 0.0f);
		cuFloatComplex training = make_cuFloatComplex(0.0f, 0.0f);

		for( unsigned int c=0; c<num_coils; c++){	
			cuFloatComplex val = csm[c*num_elements+idx];
			cuFloatComplex cval = cuConjf(val);
			norm = cuCaddf( norm, cuCmulf( val, cval ));
		    training = cuCaddf( training, cuCmulf( lowres[c*num_elements+idx], cval ));
		}

		if( cuCrealf(norm) )
			training = cuCdivf( training, norm );
	
		ctraining[idx] = training;
	}
}

void 
extract_training_data_xt( unsigned int num_coils, unsigned int num_elements, const cuFloatComplex *lowres, const cuFloatComplex *csm, cuFloatComplex *theta )
{
	dim3 blockDim(512,1,1);
	dim3 gridDim((unsigned int) ceil((double)(num_elements/blockDim.x)), 1, 1 );

	assert(num_elements>=blockDim.x);

	extract_training_data_xt_kernel<<<gridDim, blockDim>>>( num_coils, num_elements, lowres, csm, theta );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'extract_training_data_xt_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}



// Template instantiation

template bool upload_csm( mr_recon::NFFT_iteration_plan<uint2, float2, 0>*, cuFloatComplex* );
template bool upload_csm( mr_recon::NFFT_iteration_plan<uint2, float2, 1>*, cuFloatComplex* );
template bool upload_csm( mr_recon::NFFT_iteration_plan<uint2, float2, 2>*, cuFloatComplex* );

template bool upload_regularization_image( mr_recon::NFFT_iteration_plan< uint2,float2, 0>*, float*);
template bool upload_regularization_image( mr_recon::NFFT_iteration_plan< uint3,float3, 0>*, float*);
template bool upload_regularization_image( mr_recon::NFFT_iteration_plan< uint4,float4, 0>*, float*);

template bool upload_regularization_image( mr_recon::NFFT_iteration_plan< uint2,float2, 1>*, float*);
template bool upload_regularization_image( mr_recon::NFFT_iteration_plan< uint3,float3, 1>*, float*);
template bool upload_regularization_image( mr_recon::NFFT_iteration_plan< uint4,float4, 1>*, float*);

template bool upload_regularization_image( mr_recon::NFFT_iteration_plan< uint2,float2, 2>*, float*);
template bool upload_regularization_image( mr_recon::NFFT_iteration_plan< uint3,float3, 2>*, float*);
template bool upload_regularization_image( mr_recon::NFFT_iteration_plan< uint4,float4, 2>*, float*);

template bool extract_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint2, float2, 1>*, float, unsigned int, const cuFloatComplex*, float, unsigned int, unsigned int );
template bool extract_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint3, float3, 1>*, float, unsigned int, const cuFloatComplex*, float, unsigned int, unsigned int );
template bool extract_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint2, float2, 2>*, float, unsigned int, const cuFloatComplex*, float, unsigned int, unsigned int );
template bool extract_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint3, float3, 2>*, float, unsigned int, const cuFloatComplex*, float, unsigned int, unsigned int );

template bool update_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint2, float2, 1>*, const cuFloatComplex*, float, cuFloatComplex*, unsigned int, bool, unsigned int, bool );
template bool update_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint3, float3, 1>*, const cuFloatComplex*, float, cuFloatComplex*, unsigned int, bool, unsigned int, bool );
template bool update_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint2, float2, 2>*, const cuFloatComplex*, float, cuFloatComplex*, unsigned int, bool, unsigned int, bool );
template bool update_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint3, float3, 2>*, const cuFloatComplex*, float, cuFloatComplex*, unsigned int, bool, unsigned int, bool );

//template bool update_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint2, float2, 0>*, cuFloatComplex*, float2*, float*, float, cuFloatComplex*, bool, bool );
//template bool update_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint3, float3, 0>*, cuFloatComplex*, float3*, float*, float, cuFloatComplex*, bool, bool );
template bool update_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint2, float2, 0>*, cuFloatComplex*, float2*, float*, float, cuFloatComplex*, cuFloatComplex*, unsigned int, bool, bool );
template bool update_csm_and_regularization( mr_recon::NFFT_iteration_plan< uint3, float3, 0>*, cuFloatComplex*, float3*, float*, float, cuFloatComplex*, cuFloatComplex*, unsigned int, bool, bool );

template bool extract_csm_and_training_data( mr_recon::NFFT_iteration_plan< uint3, float3, 0>*, float, unsigned int, const cuFloatComplex*, float );
template bool extract_csm_and_training_data( mr_recon::NFFT_iteration_plan< uint3, float3, 1>*, float, unsigned int, const cuFloatComplex*, float );
template bool extract_csm_and_training_data( mr_recon::NFFT_iteration_plan< uint3, float3, 2>*, float, unsigned int, const cuFloatComplex*, float );

template bool extract_and_combine_current_frame( mr_recon::NFFT_iteration_plan<uint2, float2, 0>*, cuFloatComplex *, cuFloatComplex *);
template bool extract_and_combine_current_frame( mr_recon::NFFT_iteration_plan<uint2, float2, 1>*, cuFloatComplex *, cuFloatComplex *);
template bool extract_and_combine_current_frame( mr_recon::NFFT_iteration_plan<uint2, float2, 2>*, cuFloatComplex *, cuFloatComplex *);

