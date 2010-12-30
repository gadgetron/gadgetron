/*
	CUDA iterative Sense implementation.
	
	Reference:
	----------
	Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit. 
	T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
	IEEE Transactions on Medical Imaging 2009; 28(12): 1974-1985. 
*/

#include "NSense.hcu"
#include "NSense_private.hcu"
#include "preprocess_private.hcu"
#include "preprocess_radial.hcu"
#include "image_utilities.hcu"
#include "uint_util.hcu"
#include "FFT.hcu"

#include <cuComplex.h>
#include <cublas.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

using namespace mr_recon;


template< class UINTd, class FLOATd, char TYPE > NFFT_iteration_plan<UINTd, FLOATd, TYPE>*
mr_recon::preprocess_NSense_radial( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, unsigned int num_projections, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
  return preprocess_radial_NFFT<UINTd, FLOATd, TYPE>( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_samples, domain_size_coils, W, num_projections, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor );
}


/*

	Public interface

*/

// Cuda based iterative Sense: initialize

template< class UINTd, class FLOATd, char TYPE > __host__ bool
iterative_sense_initialize( NFFT_iteration_plan<UINTd,FLOATd,TYPE> *plan, unsigned int num_coils )
{
	// For now we require that we process only "full batches"
	if( num_coils % plan->domain_size_coils ){
		printf("\nNSense initialization failed! Plan's domain_size_coils must be a divisor of the number of coils.\n");
		return false;
	}

	// Initialize plan
	bool success = NFFT_initialize( plan );
	
	// If initialization failed we can't go on
	if( !success ){
		plan->initialized = false;
		printf("\nNSense initialization failed!\n");
		return false;
	}

	// Set global variables
	plan->number_of_coils = num_coils;

	return true;
}

// Cuda based iterative Sense: compute.

template< class UINTd, class FLOATd, char TYPE > __host__ bool 
iterative_sense_compute( NFFT_iteration_plan<UINTd,FLOATd,TYPE> *plan, unsigned int num_iterations, float kappa, cuFloatComplex *in_samplesDevPtr, cuFloatComplex *out_imageDevPtr )
{
	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected when starting 'iterative_sense_compute': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// We can only proceed if initialization has been done
	if( !plan->initialized ){
		printf("\n'iterative_sense_compute' failed! Call 'iterative_sense_initialize' first.\n");
		return false;
	}

	// We can only proceed if coil sensitivity maps are present
	if( !plan->CSMdevPtr ){
		printf("\n'iterative_sense_compute' failed! No coil sensitivity maps computed or uploaded.\n");
		return false;
	}

	// Ensure that kappa is not negative
	if( kappa<0 ){
		printf("\n'iterative_sense_compute': a negative kappa is not accepted.\n");
		return false;
	}

	// Number of image elements
	unsigned int num_elements = prod(plan->matrix_size);

	// Return value (we start out optimistic ;-))
	bool success = true;

	// If intensity correction image has not yet been computed then do so
	if( !plan->intensity_correction_magnitudes_image_DevPtr ){
		cudaMalloc( (void**) &plan->intensity_correction_magnitudes_image_DevPtr, num_elements*sizeof(float) );
		calculate_intensity_correction_magnitudes_image( num_elements, plan->number_of_coils, plan->CSMdevPtr, plan->intensity_correction_magnitudes_image_DevPtr );
	}

	// Allocate accumulation buffer
	cuFloatComplex *accBuffer;
	cudaMalloc( (void**) &accBuffer, num_elements*sizeof(cuFloatComplex) );

	/*
		Form initial right hand side ('accBuffer'='a')
	*/

	{
		// Clear accumulation buffer 
		clear_image( num_elements, make_cuFloatComplex(0.0f, 0.0f), accBuffer );

		// Temporary images
		cuFloatComplex *images_os, *images;
		cudaMalloc( (void**) &images, plan->domain_size_coils*num_elements*sizeof(cuFloatComplex) );
		cudaMalloc( (void**) &images_os, plan->domain_size_coils*prod(plan->matrix_size_os)*sizeof(cuFloatComplex) );

		// D
		density_compensate( plan, in_samplesDevPtr, in_samplesDevPtr );

		// Sum over coils (initialize ensures there is no remainder)
		for( unsigned int i=0; i<plan->number_of_coils/plan->domain_size_coils; i++ ){ 

			/* 
				FT1
			*/

			// Convolution

#ifdef VISUAL_TRACE
			printf("\nBefore convolution (H)"); fflush(stdout);
#endif

			success = NFFT_iteration_convolve_to_image( plan, &in_samplesDevPtr[i*plan->domain_size_coils*plan->number_of_samples], images_os );

#ifdef VISUAL_TRACE
			printf("\nAfter convolution (H)"); fflush(stdout);
#endif
			// FFT
			cuFloatComplex *tmp;
			cudaMalloc( (void**)&tmp, plan->domain_size_coils*prod(plan->matrix_size_os)*sizeof(cuFloatComplex) );

#ifdef VISUAL_TRACE
			printf("\nBefore FFT"); fflush(stdout);
#endif

			fft_shift( images_os, tmp, plan->matrix_size_os, plan->domain_size_coils);
			if( success )
				success = K2I_ALL( tmp, plan->matrix_size_os, plan->domain_size_coils, false, false );
			fft_shift( tmp, images_os, plan->matrix_size_os, plan->domain_size_coils);

			cudaFree(tmp);
			
#ifdef VISUAL_TRACE
			printf("\nAfter FFT"); fflush(stdout);
#endif

			// Deapodization
			if( success ) 
				deapodize( plan, true, images_os );

			if( success ) 
				crop_image( plan->matrix_size, plan->matrix_size_os, images, images_os, plan->domain_size_coils );

			// Sum -- S*			
			if( success )
				image_Caxpy( num_elements, true, &plan->CSMdevPtr[i*plan->domain_size_coils*num_elements], images, accBuffer, plan->domain_size_coils );
		}

		// Intensity correction
		intensity_correct_image( num_elements, prod(plan->matrix_size_os), kappa, plan->intensity_correction_magnitudes_image_DevPtr, plan->regularizationDevPtr, accBuffer );

		// We don't need the temporary images anymore
		cudaFree( images ); 
		cudaFree( images_os );
	}

	if( !success ){
		printf("\nIterative Sense failed computing 'a'!\n");
		return false;
	}

/*	
	//DEBUG
	float* tmp = (float*) calloc( 1, prod(plan->matrix_size)*sizeof(cuFloatComplex) );
	if( sizeof(UINTd)==sizeof(uint3) )
		I2K(accBuffer, uintd_to_uint4_with_ones(plan->matrix_size), 2, false);
	cudaMemcpy( tmp, accBuffer, prod(plan->matrix_size)*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("rhs.raw", "wb");
	fwrite( tmp, prod(plan->matrix_size), sizeof(cuFloatComplex), fout );
	fclose(fout);
*/	

	/*
	
		Conjugate gradient loop
		
	*/

	// CG vectors allocation
	cuFloatComplex *b, *p, *r, *q;
	cudaMalloc( (void**) &b, num_elements*sizeof(cuFloatComplex) );
	cudaMalloc( (void**) &p, num_elements*sizeof(cuFloatComplex) );
	cudaMalloc( (void**) &r, num_elements*sizeof(cuFloatComplex) );	
	cudaMalloc( (void**) &q, num_elements*sizeof(cuFloatComplex) );	

	// Copy 'a' (i.e. 'accBuffer') to 'p' and 'r'
	cudaMemcpy( p, accBuffer, num_elements*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );
	cudaMemcpy( r, accBuffer, num_elements*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

	// Clear solution (b)
	clear_image( num_elements, make_cuFloatComplex(0.0f, 0.0f), b );

	{
		// Temporary storage for the intensity corrected, and sensitivity corrected images
		cuFloatComplex *Iimage, *SIimages;

		// If we don't want to overwrite the input sample values we need a temporary array
		cuFloatComplex *cg_sample_values;

		cudaMalloc( (void**) &Iimage, num_elements*sizeof(cuFloatComplex) );
		cudaMalloc( (void**) &SIimages, plan->domain_size_coils*num_elements*sizeof(cuFloatComplex) );
		cudaMalloc( (void**) &cg_sample_values, plan->domain_size_coils*plan->number_of_samples*sizeof(cuFloatComplex) );

		// a^Ha is used for convergence checking (constant).
		const cuFloatComplex a_dot_a = cublasCdotc ( num_elements, accBuffer, 1, accBuffer, 1 );

		// r^Hr is needed throughout the CG loop (updated each iteration)
		cuFloatComplex r_dot_r = cublasCdotc ( num_elements, r, 1, r, 1 );

#ifdef PRINT_RESIDUALS
		printf("\nDelta: " ); 
#endif
		for( unsigned int iteration=0; iteration<num_iterations; iteration++ ){

			// Convergence? Report delta.
			cuFloatComplex delta = cuCdivf( r_dot_r, a_dot_a );

			// Report convergence measure
#ifdef PRINT_RESIDUALS
			printf("(%.2g) ", delta.x ); 
#endif
			/*
				NFFT iteration (E^H D E)
			*/

			// Clear accumulation buffer 
			clear_image( num_elements, make_cuFloatComplex(0.0f, 0.0f), accBuffer );

			// Intensity correction
			cudaMemcpy( Iimage, p, num_elements*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );
		    intensity_correct_image( num_elements, prod(plan->matrix_size_os), kappa, plan->intensity_correction_magnitudes_image_DevPtr, plan->regularizationDevPtr, Iimage );

#ifdef VISUAL_TRACE
		    std::cout << "BEFORE Summing over coils in cg loop" << std::endl;
#endif
			// Sum over coils (initialize guarantees that there is no remainder)
			for( unsigned int coil=0; coil<plan->number_of_coils/plan->domain_size_coils; coil++ ){

				// Sensitivity correct for each coil
				for( unsigned int _coil=0; _coil<plan->domain_size_coils; _coil++ )
					cudaMemcpy( &SIimages[_coil*num_elements], Iimage, num_elements*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );
				image_Cscal( plan->domain_size_coils*num_elements, &plan->CSMdevPtr[coil*plan->domain_size_coils*num_elements], SIimages );

#ifdef VISUAL_TRACE
		    std::cout << "BEFORE NFFT_compute" << std::endl;
#endif

				// FT2 D FT1
				if( success )
					success = NFFT_compute( plan, cg_sample_values, false, SIimages, true );

#ifdef VISUAL_TRACE
		    std::cout << "AFTER NFFT_compute" << std::endl;
#endif
				// Sum -- S*			
				if( success )
					image_Caxpy( num_elements, true, &plan->CSMdevPtr[coil*plan->domain_size_coils*num_elements], SIimages, accBuffer, plan->domain_size_coils );
			}			
			
#ifdef VISUAL_TRACE
		    std::cout << "AFTER Summing over coils in cg loop" << std::endl;
#endif
			// Regularization
			if( plan->regularizationDevPtr )
				image_Cssaxpy( num_elements, kappa, plan->regularizationDevPtr, Iimage, accBuffer );

			// Intensity correction
			intensity_correct_image( num_elements, prod(plan->matrix_size_os), kappa, plan->intensity_correction_magnitudes_image_DevPtr, plan->regularizationDevPtr, accBuffer );

			// Copy image to 'q'
			cudaMemcpy( q, accBuffer, num_elements*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

			// Scaling dot products
			cuFloatComplex p_dot_q = cublasCdotc ( num_elements, p, 1, q, 1 );
			cuFloatComplex scale = cuCdivf( r_dot_r, p_dot_q );

			// Update 'b'
			cublasCaxpy ( num_elements, scale, p, 1, b, 1 );

			// update 'r'
			scale.x *= -1.0f; scale.y *= -1.0f;
			cublasCaxpy( num_elements, scale, q, 1, r, 1 );

			// Temporarily save old r_dot_r and update it
			cuFloatComplex r_dot_r_prev = r_dot_r;
			r_dot_r = cublasCdotc ( num_elements, r, 1, r, 1 );
			scale = cuCdivf( r_dot_r, r_dot_r_prev );

			// Update 'p'
			cublasCscal( num_elements, scale, p, 1 );
			cublasCaxpy( num_elements, make_cuFloatComplex(1.0f, 0.0f), r, 1, p, 1 );

/*
	//DEBUG
	float* tmp = (float*) calloc( 1, prod(plan->matrix_size)*sizeof(cuFloatComplex) );
	if( sizeof(UINTd)==sizeof(uint3) )
		I2K(accBuffer, uintd_to_uint4_with_ones(plan->matrix_size), 2, false);
	cudaMemcpy( tmp, accBuffer, prod(plan->matrix_size)*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("iteration.raw", "wb");
	fwrite( tmp, prod(plan->matrix_size), sizeof(cuFloatComplex), fout );
	fclose(fout);
*/
		}

		// We don't need the temporary images anymore
		cudaFree( Iimage ); 
		cudaFree( SIimages ); 
		cudaFree( cg_sample_values );
	}

	// Final intensity correction
	intensity_correct_image( num_elements, prod(plan->matrix_size_os), kappa, plan->intensity_correction_magnitudes_image_DevPtr, plan->regularizationDevPtr, b );
	cudaMemcpy( out_imageDevPtr, b, num_elements*sizeof(cuFloatComplex), cudaMemcpyDeviceToDevice );

/*
	//DEBUG
	printf("\nWriting result to disk"); fflush(stdout);
	float* tmp = (float*) calloc( 1, prod(plan->matrix_size)*sizeof(cuFloatComplex) );
	I2K(b, uintd_to_uint4_with_ones(plan->matrix_size), 2, false);
	float *_devPtr;
	cudaMalloc( (void**) &_devPtr, prod(plan->matrix_size)*sizeof(float) );
	image_modulus( b, _devPtr, prod(plan->matrix_size), true );
	cudaMemcpy( tmp, _devPtr, prod(plan->matrix_size)*sizeof(float), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("result.raw", "wb");
	fwrite( tmp, prod(plan->matrix_size), sizeof(float), fout );
	fclose(fout);
*/

	// Cleanup
	cudaFree(b); 
	cudaFree(p); 
	cudaFree(r); 
	cudaFree(q);
	cudaFree( accBuffer );

	return success;
}

/*

	Private interface

*/


// Calculate intensity correction image
__host__ void
calculate_intensity_correction_magnitudes_image( unsigned int num_elements, unsigned int num_coils, cuFloatComplex *coilMapsDevPtr, float *intensityImageDevPtr )
{

	// Find dimensions of grid/blocks (this can be very small problems).

	dim3 dimBlock( 512, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	calculate_intensity_correction_magnitudes_image_kernel<<< dimGrid, dimBlock >>>( num_elements, num_coils, coilMapsDevPtr, intensityImageDevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'calculate_intensity_correction_magnitudes_image': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

// Image "Caxpy": y = ax+y where _all_ of a,x, and y are vectors of complex float. 'y' is a single image (accumulation buffer).
__host__ void 
image_Caxpy( unsigned int num_elements, bool complexConjugate_a, cuFloatComplex *a, cuFloatComplex *x, cuFloatComplex *y, unsigned int number_of_images )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/deviceProp.maxThreadsPerBlock), 1, 1 );

	// Invoke kernel
	image_Caxpy_kernel<<< dimGrid, dimBlock >>>( num_elements, complexConjugate_a, a, x, y, number_of_images );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_Caxpy': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

// Image "Caxpy": y = ax+y where a is a scalar vector and x,y are vectors of complex float.
__host__ 
void image_Cssaxpy( unsigned int num_elements, float _a, float *a, cuFloatComplex *x, cuFloatComplex *y )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/deviceProp.maxThreadsPerBlock), 1, 1 );

	// Invoke kernel
	image_Cssaxpy_kernel<<< dimGrid, dimBlock >>>( num_elements, _a, a, x, y );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_Caxpy': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

// Intensity correct
__host__ void 
intensity_correct_image( unsigned int num_elements, unsigned int num_elements_os, float kappa, float *intensityMagnitudesImage, float *regularizationImage, cuFloatComplex *image )
{
	if( !intensityMagnitudesImage || !image ){
		printf("\nERROR: NULL pointer passed to 'intensity_correct_image'. Quitting.\n");
		exit(1);
	}

	// Find dimensions of grid/blocks (small problem).

	dim3 dimBlock( 512, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	intensity_correct_image_kernel<<< dimGrid, dimBlock >>>( num_elements, num_elements_os, kappa, intensityMagnitudesImage, regularizationImage, image );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'intensity_correct_image': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

// Multiply two complex images
__host__ void 
image_Cscal( unsigned int num_elements, cuFloatComplex *a, cuFloatComplex *x )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/deviceProp.maxThreadsPerBlock), 1, 1 );

	// Invoke kernel
	image_Cscal_kernel<<< dimGrid, dimBlock >>>( num_elements, a, x );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_Cscal': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

// Multiply float image with complex image
__host__ void 
image_Csscal( unsigned int num_elements, float *a, cuFloatComplex *x )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/deviceProp.maxThreadsPerBlock), 1, 1 );

	// Invoke kernel
	image_Csscal_kernel<<< dimGrid, dimBlock >>>( num_elements, a, x );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_Cscal': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


/*

	Kernels

*/

// Calculate intensity correction image kernel
__global__ void
calculate_intensity_correction_magnitudes_image_kernel( unsigned int num_elements, unsigned int num_coils, cuFloatComplex *coilMapsDevPtr, float *intensityImageDevPtr )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){
		float val = 0.0f;
		for( unsigned int i=0; i<num_coils; i++ ){
			cuFloatComplex z = coilMapsDevPtr[i*num_elements+idx];
			val += (z.x*z.x+z.y*z.y);
		}
		intensityImageDevPtr[idx] = val;
	}
}

// Intensity correct image (preconditioning)
__global__ void 
intensity_correct_image_kernel( unsigned int num_elements, unsigned int num_elements_os, float kappa, float *intensityMagnitudesImage, float *regularizationImage, cuFloatComplex *image )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){

		float val = intensityMagnitudesImage[idx];

		val = val*num_elements_os + ((regularizationImage) ? kappa*regularizationImage[idx] : 0.0f);

		float sqt_val = sqrtf(val);

		if( sqt_val<0.000001f )
			sqt_val = 0.000001f;
		
		cuFloatComplex pixel = image[idx];
		float scale = 1.0f/sqt_val;
		image[idx] = make_cuFloatComplex( scale*cuCrealf(pixel), scale*cuCimagf(pixel) );
	}
}
	
// Image "Caxpy": y = ax+y where _all_ of a,x, and y are vectors of complex float. 'y' is a single image (accumulation buffer).
__global__ void 
image_Caxpy_kernel( unsigned int num_elements, bool complexConjugate_a, cuFloatComplex *a, cuFloatComplex *x, cuFloatComplex *y, unsigned int number_of_images )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){
		cuFloatComplex acc = make_cuFloatComplex( 0.0f, 0.0f );
		for( unsigned int i=0; i<number_of_images; i++ ){
			cuFloatComplex _a = a[i*num_elements+idx];
			if( complexConjugate_a ) 
				_a = cuConjf(_a);
			acc = cuCaddf( acc, cuCmulf( _a, x[i*num_elements+idx] ));
		}
		y[idx] = cuCaddf( y[idx], acc );
	}
}

// Image "Cssaxpy": y = ax+y where x and y are vectors of complex float, a is float vector.
__global__ void 
image_Cssaxpy_kernel( unsigned int num_elements, float s, float *a, cuFloatComplex *x, cuFloatComplex *y )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){
		cuFloatComplex _x = x[idx];
		float _a = s*a[idx];
		y[idx] = cuCaddf( y[idx],  make_cuFloatComplex( _a*_x.x, _a*_x.y ));
	}
}

// Multiply two complex images
__global__ void 
image_Cscal_kernel( unsigned int num_elements, cuFloatComplex *a, cuFloatComplex *x )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){
		x[idx] = cuCmulf( a[idx], x[idx] );
	}
}

// Multiply float image with complex image
__global__ void 
image_Csscal_kernel( unsigned int num_elements, float *a, cuFloatComplex *x )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){
		float _a = a[idx];
		cuFloatComplex _x = x[idx];
		x[idx] = make_cuFloatComplex( _a*_x.x, _a*_x.y );
	}
}



// Instatiation

template mr_recon::NFFT_iteration_plan< uint2, float2, 1>* mr_recon::preprocess_NSense_radial< uint2,float2, 1>( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_iteration_plan< uint3, float3, 1>* mr_recon::preprocess_NSense_radial< uint3,float3, 1>( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template mr_recon::NFFT_iteration_plan< uint2, float2, 2>* mr_recon::preprocess_NSense_radial< uint2,float2, 2>( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );
template mr_recon::NFFT_iteration_plan< uint3, float3, 2>* mr_recon::preprocess_NSense_radial< uint3,float3, 2>( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float );

template  bool iterative_sense_initialize< uint2,float2, 0>( NFFT_iteration_plan< uint2,float2, 0>*, unsigned int );
template  bool iterative_sense_initialize< uint3,float3, 0>( NFFT_iteration_plan< uint3,float3, 0>*, unsigned int );
template  bool iterative_sense_initialize< uint4,float4, 0>( NFFT_iteration_plan< uint4,float4, 0>*, unsigned int );

template  bool iterative_sense_initialize< uint2,float2, 1>( NFFT_iteration_plan< uint2,float2, 1>*, unsigned int );
template  bool iterative_sense_initialize< uint3,float3, 1>( NFFT_iteration_plan< uint3,float3, 1>*, unsigned int );
template  bool iterative_sense_initialize< uint4,float4, 1>( NFFT_iteration_plan< uint4,float4, 1>*, unsigned int );

template  bool iterative_sense_initialize< uint2,float2, 2>( NFFT_iteration_plan< uint2,float2, 2>*, unsigned int );
template  bool iterative_sense_initialize< uint3,float3, 2>( NFFT_iteration_plan< uint3,float3, 2>*, unsigned int );
template  bool iterative_sense_initialize< uint4,float4, 2>( NFFT_iteration_plan< uint4,float4, 2>*, unsigned int );

template  bool iterative_sense_compute< uint2,float2, 0>( NFFT_iteration_plan< uint2,float2, 0>*, unsigned int, float, cuFloatComplex*, cuFloatComplex* );
template  bool iterative_sense_compute< uint3,float3, 0>( NFFT_iteration_plan< uint3,float3, 0>*, unsigned int, float, cuFloatComplex*, cuFloatComplex* );
template  bool iterative_sense_compute< uint4,float4, 0>( NFFT_iteration_plan< uint4,float4, 0>*, unsigned int, float, cuFloatComplex*, cuFloatComplex* );

template  bool iterative_sense_compute< uint2,float2, 1>( NFFT_iteration_plan< uint2,float2, 1>*, unsigned int, float, cuFloatComplex*, cuFloatComplex* );
template  bool iterative_sense_compute< uint3,float3, 1>( NFFT_iteration_plan< uint3,float3, 1>*, unsigned int, float, cuFloatComplex*, cuFloatComplex* );
template  bool iterative_sense_compute< uint4,float4, 1>( NFFT_iteration_plan< uint4,float4, 1>*, unsigned int, float, cuFloatComplex*, cuFloatComplex* );

template  bool iterative_sense_compute< uint2,float2, 2>( NFFT_iteration_plan< uint2,float2, 2>*, unsigned int, float, cuFloatComplex*, cuFloatComplex* );
template  bool iterative_sense_compute< uint3,float3, 2>( NFFT_iteration_plan< uint3,float3, 2>*, unsigned int, float, cuFloatComplex*, cuFloatComplex* );
template  bool iterative_sense_compute< uint4,float4, 2>( NFFT_iteration_plan< uint4,float4, 2>*, unsigned int, float, cuFloatComplex*, cuFloatComplex* );
