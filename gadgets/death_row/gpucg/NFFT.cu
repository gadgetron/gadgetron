/*

	CUDA implementation of the NFFT and NFFT_H algorithms (gridding and inverse gridding).


	References:
	-----------

	Accelerating the Non-equispaced Fast Fourier Transform on Commodity Graphics Hardware.
	T.S. Sørensen, T. Schaeffter, K.Ø. Noe, M.S. Hansen. 
	IEEE Transactions on Medical Imaging 2008; 27(4):538-547.

	Real-time Reconstruction of Sensitivity Encoded Radial Magnetic Resonance Imaging Using a Graphics Processing Unit.
	T.S. Sørensen, D. Atkinson, T. Schaeffter, M.S. Hansen.
	IEEE Transactions on Medical Imaging 2009; 28(12): 1974-1985. 
*/


// Includes - our own code
#include "NFFT.hcu"
#include "NFFT_private.hcu"
#include "preprocess.hpp"
#include "preprocess_private.hcu"
#include "image_utilities.hcu"
#include "int_util.hcu"
#include "uint_util.hcu"
#include "uint_util_device.hcu"
#include "float_util.hcu"
#include "float_util_device.hcu"
#include "float_util_host.hcu"
#include "FFT.hcu"

// Includes - CUDA
#include <math_constants.h>
#include <cuComplex.h>
#include <cufft.h>
#include <cublas.h>

// Includes - stdlibs
#include <stdio.h>
#include <assert.h>
#include <limits.h>

// This is the device number we will use 
unsigned int _convolution_device = 0;

/*
       Some defines used to configure our kernels
*/

// Important to ensure stability!!!
#define BUFFER_SAFETY_CHECK

// It seems from the occupancy calculator that the G80 allocates shared memory in 512-byte chunks
#define MIN_SHARED_MEM_PER_BLOCK 512


// Block configuration for kernel
// ------------------------------
// gen = generic, sgen = semigeneric, par = parallel, seq = sequential
// x,y -> w threads means: shared memory = x, registers = y, max threads = w
// This information is obtained from the .cubin file (x,y) and the Occupancy Calculator (w).


// gen: 80,32 -> 256 threads
// par: 80,29 -> 256 threads	(24->320 with --maxrregcount)
// seq: 80,28 -> 256 threads

#define NFFT_SHARED_MEMORY_OVERHEAD_PER_2D_KERNEL_BLOCK	80
#define NFFT_THREADS_PER_2D_KERNEL_GEN						256
#define NFFT_THREADS_PER_2D_KERNEL_PAR						256
#define NFFT_THREADS_PER_2D_KERNEL_SEQ						256

// gen: 80,36 -> 192 threads
// par: 80,32 -> 256 threads
// seq: 80,34 -> 192 threads

#define NFFT_SHARED_MEMORY_OVERHEAD_PER_3D_KERNEL_BLOCK	80
#define NFFT_THREADS_PER_3D_KERNEL_GEN						192
#define NFFT_THREADS_PER_3D_KERNEL_PAR_G80					192
#define NFFT_THREADS_PER_3D_KERNEL_PAR_G200					384
#define NFFT_THREADS_PER_3D_KERNEL_SEQ						192

// gen: 80,40 -> 192 threads
// par: 80,37 -> 192 threads
// seq: 80,38 -> 192 threads

#define NFFT_SHARED_MEMORY_OVERHEAD_PER_4D_KERNEL_BLOCK	80
#define NFFT_THREADS_PER_4D_KERNEL_GEN						192
#define NFFT_THREADS_PER_4D_KERNEL_PAR						192
#define NFFT_THREADS_PER_4D_KERNEL_SEQ						192

// gen: 96,29 -> 256 threads
// sgen: 96,25 -> 320 threads
// par: 96,23 -> 320 threads
// seq: 96,26 -> 256 threads (24->320 with --maxrregcount)

#define NFFT_H_SHARED_MEMORY_OVERHEAD_PER_2D_KERNEL_BLOCK	96
#define NFFT_H_THREADS_PER_2D_KERNEL_GEN					256
#define NFFT_H_THREADS_PER_2D_KERNEL_SEMI_GEN				320
#define NFFT_H_THREADS_PER_2D_KERNEL_PAR					512
#define NFFT_H_THREADS_PER_2D_KERNEL_SEQ					256

// gen: 112,35 -> 192 threads
// sgen: 112,27 -> 256 threads (24->320 with --maxrregcount)
// par: 96,24 -> 320 threads
// seq: 112,32 -> 192 threads

#define NFFT_H_SHARED_MEMORY_OVERHEAD_PER_3D_KERNEL_BLOCK	112
#define NFFT_H_THREADS_PER_3D_KERNEL_GEN					192
#define NFFT_H_THREADS_PER_3D_KERNEL_SEMI_GEN				256
//#define NFFT_H_THREADS_PER_3D_KERNEL_PAR_G80				320
//#define NFFT_H_THREADS_PER_3D_KERNEL_PAR					256
#define NFFT_H_THREADS_PER_3D_KERNEL_PAR_G80				256
#define NFFT_H_THREADS_PER_3D_KERNEL_PAR_G200				512
#define NFFT_H_THREADS_PER_3D_KERNEL_SEQ					192

// gen: 144,42 -> 192 threads
// sgen: 144,29 -> 256 threads
// par: 128,25 -> 320 threads
// seq: 128,38 -> 192 threads

#define NFFT_H_SHARED_MEMORY_OVERHEAD_PER_4D_KERNEL_BLOCK	144
#define NFFT_H_THREADS_PER_4D_KERNEL_GEN					192
#define NFFT_H_THREADS_PER_4D_KERNEL_SEMI_GEN				256
#define NFFT_H_THREADS_PER_4D_KERNEL_PAR					320
#define NFFT_H_THREADS_PER_4D_KERNEL_SEQ					192


/*
       Declare the textures required to implement convolution. 
       Use linear global memory.
*/

// NFFT_iteration (to compute E^H*E)
texture<float2, 1> tex_sample_positions_f2_it;		// For 2D trajectories
texture<float4, 1> tex_sample_positions_f4_it;		// For 4D trajectories (used also for 3D traj.)
texture<cuFloatComplex, 1> tex_grid_values_it;		// Complex sample values
texture<cuFloatComplex, 1> tex_sample_values_it;	// Complex sample values
texture<uint2, 1> tex_stripsMap_NFFT_it;		// Strips map (maps domain indices to strips)
texture<uint2, 1> tex_stripsDir_ui2_NFFT_it;		// Strips encoding direction for 2D grids
texture<uint4, 1> tex_stripsDir_ui4_NFFT_it;		// Strips encoding direction for 4D grids (used also for 3D grids)
texture<uint2, 1> tex_stripOrigins_ui2_NFFT_it;		// Strips origins for 2D grids
texture<uint4, 1> tex_stripOrigins_ui4_NFFT_it;		// Strips origins for 4D grids (used also for 3D grids)
texture<unsigned int, 1> tex_stripLengths_NFFT_it;	// Strip lenghts
texture<unsigned int, 1> tex_domainsMap_NFFT_H_it;	// Domains map (maps thread indices to domains)
texture<uint2, 1> tex_stripsMap_NFFT_H_it;		// Strips map (maps domain indices to strips)
texture<uint2, 1> tex_strips_NFFT_H_it;			// Strips

// NFFT_H (gridding)
texture<float2, 1> tex_sample_positions_f2_NFFT_H;	// For 2D trajectories
texture<float4, 1> tex_sample_positions_f4_NFFT_H;	// For 4D trajectories (used also for 3D traj.)
texture<cuFloatComplex, 1> tex_sample_values_NFFT_H;	// Complex sample values
texture<unsigned int, 1> tex_domainsMap_NFFT_H;		// Domains map (maps thread indices to domains)
texture<uint2, 1> tex_stripsMap_NFFT_H;			// Strips map (maps domain indices to strips)
texture<uint2, 1> tex_strips_NFFT_H;			// Strips

// WARNING: Compiler issue
// Some of the textures below will be numbered >32 i.e. ILLEGALLY!!!
// We don't need them for now so it doesn't matter (yet)...

// NFFT (inverse gridding)
texture<float2, 1> tex_sample_positions_f2_NFFT;	// For 2D trajectories
texture<float4, 1> tex_sample_positions_f4_NFFT;	// For 4D trajectories (used also for 3D traj.)
texture<cuFloatComplex, 1> tex_grid_values_NFFT;	// Complex sample values
texture<uint2, 1> tex_stripsMap_NFFT;			// Strips map (maps domain indices to strips)
texture<uint2, 1> tex_stripsDir_ui2_NFFT;		// Strips encoding direction for 2D grids
texture<uint4, 1> tex_stripsDir_ui4_NFFT;		// Strips encoding direction for 4D grids (used also for 3D grids)
texture<uint2, 1> tex_stripOrigins_ui2_NFFT;		// Strips origins for 2D grids
texture<uint4, 1> tex_stripOrigins_ui4_NFFT;		// Strips origins for 4D grids (used also for 3D grids)
texture<unsigned int, 1> tex_stripLengths_NFFT;		// Strip lenghts

/*
	The textures above are accessed with these macros.
*/

// NFFT
#define TEX_SAMPLE_POSITIONS_F2_NFFT(b) (b) ? tex_sample_positions_f2_NFFT : tex_sample_positions_f2_it
#define TEX_SAMPLE_POSITIONS_F4_NFFT(b) (b) ? tex_sample_positions_f4_NFFT : tex_sample_positions_f4_it
#define TEX_GRID_VALUES_NFFT(b) (b) ? tex_grid_values_NFFT : tex_grid_values_it
#define TEX_STRIPS_MAP_NFFT(b) (b) ? tex_stripsMap_NFFT : tex_stripsMap_NFFT_it
#define TEX_STRIPS_DIR_UI2_NFFT(b) (b) ? tex_stripsDir_ui2_NFFT : tex_stripsDir_ui2_NFFT_it
#define TEX_STRIPS_DIR_UI4_NFFT(b) (b) ? tex_stripsDir_ui4_NFFT : tex_stripsDir_ui4_NFFT_it
#define TEX_STRIPS_ORIGINS_UI2_NFFT(b) (b) ? tex_stripOrigins_ui2_NFFT : tex_stripOrigins_ui2_NFFT_it
#define TEX_STRIPS_ORIGINS_UI4_NFFT(b) (b) ? tex_stripOrigins_ui4_NFFT : tex_stripOrigins_ui4_NFFT_it
#define TEX_STRIPS_LENGTHS_NFFT(b) (b) ? tex_stripLengths_NFFT : tex_stripLengths_NFFT_it

// NFFT_H
#define TEX_SAMPLE_POSITIONS_F2_NFFT_H(b) (b) ? tex_sample_positions_f2_NFFT_H : tex_sample_positions_f2_it
#define TEX_SAMPLE_POSITIONS_F4_NFFT_H(b) (b) ? tex_sample_positions_f4_NFFT_H : tex_sample_positions_f4_it
#define TEX_SAMPLE_VALUES_NFFT_H(b) (b) ? tex_sample_values_NFFT_H : tex_sample_values_it
#define TEX_DOMAINS_MAP_NFFT_H(b) (b) ? tex_domainsMap_NFFT_H : tex_domainsMap_NFFT_H_it
#define TEX_STRIPS_MAP_NFFT_H(b) (b) ? tex_stripsMap_NFFT_H : tex_stripsMap_NFFT_H_it
#define TEX_STRIPS_NFFT_H(b) (b) ? tex_strips_NFFT_H : tex_strips_NFFT_H_it

/*
        Reference to shared memory
*/

extern __shared__ /*volatile*/ float shared_mem[];


/*
	Variables specific for the "online radial" types
*/

__constant__ float __num_samples_per_projection;
__constant__ float __num_samples_per_frame;
__constant__ float __num_projections_per_frame;
__constant__ float __one_over_num_samples_per_projection;
__constant__ float __one_over_num_samples_per_frame;
__constant__ float __one_over_num_projections_per_frame;
__constant__ float __angular_offset;
__constant__ float __gc_factor;
__constant__ float __frames_per_rotation_cycle;
__constant__ float __rotation_gap;
__constant__ float __interframe_rotation;
__constant__ float __one_over_radial_oversampling;

#ifdef SELF_SCHEDULING

// Warp size
#define WARP_SIZE 32

// How many warps to put in each block (blockDim.x==WARP_SIZE)
#define BLOCKDIM_Y 6

__constant__ unsigned int globalPoolThreadCount;
__device__ unsigned int globalPoolNextThread;

#endif


template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ void _nfft_NFFT_set_constants( PLAN<UINTd, FLOATd, TYPE> *plan )
{
	// Utility to set constant memory space

	float _one_over_radial_oversampling = (float)plan->matrix_size.x/(float)plan->samples_per_projection;
	cudaMemcpyToSymbol( __one_over_radial_oversampling, &_one_over_radial_oversampling, sizeof(float) );

	float _num_samples_per_projection = (float)plan->samples_per_projection;
	cudaMemcpyToSymbol( __num_samples_per_projection, &_num_samples_per_projection, sizeof(float) );

	float _num_samples_per_frame = (float)(plan->samples_per_projection*plan->projections_per_frame);
	cudaMemcpyToSymbol( __num_samples_per_frame, &_num_samples_per_frame, sizeof(float) );

	float _one_over_num_samples_per_frame = 1.0f/(float)(plan->samples_per_projection*plan->projections_per_frame);
	cudaMemcpyToSymbol( __one_over_num_samples_per_frame, &_one_over_num_samples_per_frame, sizeof(float) );

	float _one_over_num_samples_per_projection = 1.0f/(float)plan->samples_per_projection;
	cudaMemcpyToSymbol( __one_over_num_samples_per_projection, &_one_over_num_samples_per_projection, sizeof(float) );

	float _one_over_num_projections_per_frame = 1.0f/(float)plan->projections_per_frame;
	cudaMemcpyToSymbol( __one_over_num_projections_per_frame, &_one_over_num_projections_per_frame, sizeof(float) );

	float _num_projections_per_frame = (float)(plan->projections_per_frame);
	cudaMemcpyToSymbol( __num_projections_per_frame, &_num_projections_per_frame, sizeof(float) );

	switch(TYPE){
	
		case 1:
			{
				cudaMemcpyToSymbol( __gc_factor, &plan->gc_factor, sizeof(float) );

				float _angular_offset = (float)plan->angular_offset;
				cudaMemcpyToSymbol( __angular_offset, &_angular_offset, sizeof(float) );
			}
			break;

		case 2:
			{
				cudaMemcpyToSymbol( __interframe_rotation, &plan->interframe_rotation, sizeof(float) );

				float _frames_per_rotation_cycle = (float)plan->frames_per_rotation;
				cudaMemcpyToSymbol( __frames_per_rotation_cycle, &_frames_per_rotation_cycle, sizeof(float) );

				float _rotation_gap = (float)((plan->frames_per_rotation-1)/2);
				cudaMemcpyToSymbol( __rotation_gap, &_rotation_gap, sizeof(float) );

				float _angular_offset = (float)plan->angular_offset;
				cudaMemcpyToSymbol( __angular_offset, &_angular_offset, sizeof(float) );
			}
			break;

		default:
			printf("\nInternal ERROR!: 'NFFT_set_constants': switch. Quitting.\n" ); fflush(stdout); 
			exit(1);
			break;
	}	

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nInternal ERROR!: 'NFFT_set_constants': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}
}

template< char TYPE > 
__host__ void _nfft_NFFT_set_constants( unsigned int matrix_size, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor )
{
	// Utility to set constant memory space

	float _one_over_radial_oversampling = (float)matrix_size/(float)samples_per_projection;
	cudaMemcpyToSymbol( __one_over_radial_oversampling, &_one_over_radial_oversampling, sizeof(float) );

	float _num_samples_per_projection = (float)samples_per_projection;
	cudaMemcpyToSymbol( __num_samples_per_projection, &_num_samples_per_projection, sizeof(float) );

	float _num_samples_per_frame = (float)(samples_per_projection*projections_per_frame);
	cudaMemcpyToSymbol( __num_samples_per_frame, &_num_samples_per_frame, sizeof(float) );

	float _one_over_num_samples_per_frame = 1.0f/(float)_num_samples_per_frame;
	cudaMemcpyToSymbol( __one_over_num_samples_per_frame, &_one_over_num_samples_per_frame, sizeof(float) );

	float _one_over_num_samples_per_projection = 1.0f/(float)samples_per_projection;
	cudaMemcpyToSymbol( __one_over_num_samples_per_projection, &_one_over_num_samples_per_projection, sizeof(float) );

	float _one_over_num_projections_per_frame = 1.0f/(float)projections_per_frame;
	cudaMemcpyToSymbol( __one_over_num_projections_per_frame, &_one_over_num_projections_per_frame, sizeof(float) );

	float _num_projections_per_frame = (float)(projections_per_frame);
	cudaMemcpyToSymbol( __num_projections_per_frame, &_num_projections_per_frame, sizeof(float) );

	switch(TYPE){
	
		case 1:
			{
				cudaMemcpyToSymbol( __gc_factor, &gc_factor, sizeof(float) );

				float _angular_offset = (float)angular_offset;
				cudaMemcpyToSymbol( __angular_offset, &_angular_offset, sizeof(float) );
			}
			break;

		case 2:
			{
				float _interframe_rotation = CUDART_PI_F/(float)(frames_per_rotation*projections_per_frame);
				cudaMemcpyToSymbol( __interframe_rotation, &_interframe_rotation, sizeof(float) );

				float _frames_per_rotation_cycle = (float)frames_per_rotation;
				cudaMemcpyToSymbol( __frames_per_rotation_cycle, &_frames_per_rotation_cycle, sizeof(float) );

				float _rotation_gap = (float)((frames_per_rotation-1)/2);
				cudaMemcpyToSymbol( __rotation_gap, &_rotation_gap, sizeof(float) );

				float _angular_offset = (float)angular_offset;
				cudaMemcpyToSymbol( __angular_offset, &_angular_offset, sizeof(float) );
			}
			break;

		default:
			printf("\nInternal ERROR!: 'NFFT_set_constants': switch. Quitting.\n" ); fflush(stdout); 
			exit(1);
			break;
	}	

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nInternal ERROR!: 'NFFT_set_constants': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}
}


// Utility function to compute radial trajectories. Local!: depends on constant memory set initially.
template <class UINTd, class FLOATd, char TYPE> __inline__ __device__ FLOATd
cuda_compute_radial_sample_position( unsigned int sampleIdx, float alpha, UINTd bias, UINTd bias_os );

/*
        Kernel control parameters.
	TODO: Read from lookup table based on alpha and W.
       ------
*/

// Kaiser-Beseel control parameter (set in PLAN::initialize)
static float beta = 0.0f;


/*
	NFFT  - host interface implementation
*/

// Initialize NFFT

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ bool NFFT_initialize( PLAN<UINTd, FLOATd, TYPE> *plan )
{
	return plan->initialize();
}

// Compute entire NFFT

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ bool NFFT_compute( PLAN<UINTd, FLOATd, TYPE> *plan, cuFloatComplex *samplesDevPtr, bool oversampled_image, cuFloatComplex *imageDevPtr, bool densityCompensate )
{
	return plan->compute( samplesDevPtr, imageDevPtr, oversampled_image, densityCompensate );
}

// Cleanup after NFFT

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ bool NFFT_cleanup( PLAN<UINTd, FLOATd, TYPE> **plan )
{
	bool success = (*plan)->cleanup();
	*plan = 0x0;
	return success;
}


/*
	NFFT  - host building blocks interface implementation
*/


// NFFT convolution

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ bool NFFT_convolve( PLAN<UINTd, FLOATd, TYPE> *plan, cuFloatComplex *samplesDevPtr, cuFloatComplex *imageDevPtr )
{
	return plan->convolve( samplesDevPtr, imageDevPtr );
}

template< class UINTd, class FLOATd, char TYPE > 
__host__ bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE> *plan, cuFloatComplex *samplesDevPtr, cuFloatComplex *imageDevPtr )
{
	return NFFT_convolve_to_samples<UINTd, FLOATd, TYPE, false>( plan, &tex_grid_values_it, samplesDevPtr, imageDevPtr );
}

template< class UINTd, class FLOATd, char TYPE > 
__host__ bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE> *plan, cuFloatComplex *samplesDevPtr, cuFloatComplex *imageDevPtr, bool accumulate )
{
  return NFFT_convolve_to_image<UINTd, FLOATd, TYPE, false>( plan, &tex_sample_values_it, samplesDevPtr, imageDevPtr, accumulate );
}


// Weights normalization

__host__ bool
cuda_normalize_weights( unsigned int num_elements, unsigned int fixed_dims_multiplum, float area, float *weightsDevPtr )
{
	// Density compensation weights must be scaled to the area of the covered k-space.

	if( num_elements%fixed_dims_multiplum ){
		printf("\nError normalizing density compensation weights: total number of elements is not a multiple of the fixed dimensions sizes!\n");
		return false;
	}

	unsigned int elements_per_sum = num_elements/fixed_dims_multiplum;

	for( unsigned i=0; i<fixed_dims_multiplum; i++ ){

		float sum = cublasSasum( elements_per_sum, &weightsDevPtr[i*elements_per_sum], 1 );
		float scale = area / sum;

		assert( sum > 0 );

//		printf("\nDensity compensation sum: %f, scale: %f, elements_per_sum: %d", sum, scale, elements_per_sum );

		//Scale the array
		cublasSscal( elements_per_sum, scale, &weightsDevPtr[i*elements_per_sum], 1 );
	}

	//DEBUG
//	float* tmp = (float*) calloc( 1, num_elements*sizeof(cuFloatComplex) );
//	cudaMemcpy( tmp, accBuffer, prod(plan->matrix_size)*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
//	FILE *fout = fopen("weights.raw", "wb");
//	fwrite( tmp, prod(plan->matrix_size), sizeof(cuFloatComplex), fout );
//	fclose(fout);

	return true;
}


// Density compensation

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > __host__ bool 
density_compensate( PLAN<UINTd, FLOATd, TYPE> *plan, cuFloatComplex *out_samples_devPtr, const cuFloatComplex *in_samples_devPtr )
{
	if( !plan->weights_DevPtr ){
		printf("\nWarning: 'density_compensate' : density compensation weights have not been computed/uploaded. Returning.");
		return false;
	}

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );
	
	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)plan->number_of_samples/deviceProp.maxThreadsPerBlock), 1, 1 );

	assert( plan->number_of_samples >= (unsigned int)deviceProp.maxThreadsPerBlock );

	// Invoke kernel
	density_compensate_kernel<<< dimGrid, dimBlock >>> ( plan->number_of_samples, out_samples_devPtr, in_samples_devPtr, plan->weights_DevPtr, plan->number_of_coils );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'density_compensate_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}

	return true;
}


// Deapodization 

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ void deapodize( PLAN<UINTd, FLOATd, TYPE> *plan, bool oversampled_image, cuFloatComplex *imageDevPtr, bool ignore_fixed_dims )
{
	// Get device properties
	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	// The different plan types are identified based on their size. 
	assert( sizeof(mr_recon::NFFT_plan<UINTd,FLOATd,TYPE>) != sizeof(mr_recon::NFFT_H_plan<UINTd,FLOATd,TYPE>) );
	assert( sizeof(mr_recon::NFFT_plan<UINTd,FLOATd,TYPE>) != sizeof(mr_recon::NFFT_iteration_plan<UINTd,FLOATd,TYPE>) );
	assert( sizeof(mr_recon::NFFT_H_plan<UINTd,FLOATd,TYPE>) != sizeof(mr_recon::NFFT_iteration_plan<UINTd,FLOATd,TYPE>) );

	// Invoke kernel
	if( !oversampled_image ){
	  
		// Simple case: Each image matches the deapodization filter in size
	
		// If the last dimensions are fixed then the dimensionality of deapodization filter is reduced
		unsigned int fixed_dims_multiplier = 1, num_elements = 1;
		for( unsigned int _d=0; _d<plan->d; _d++ ){
			if( ((unsigned int*)&plan->fixed_dims)[_d] )
				fixed_dims_multiplier *= ((unsigned int*)&plan->matrix_size)[_d];
			else 
				num_elements *= ((unsigned int*)&plan->matrix_size)[_d];
		}

		if( ignore_fixed_dims )
			fixed_dims_multiplier = 1;

		dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
		dim3 dimGrid( (unsigned int) ceil((double)num_elements/deviceProp.maxThreadsPerBlock), 1, 1 );

		switch( sizeof(*plan))
		{
		case sizeof( mr_recon::NFFT_plan<UINTd,FLOATd,TYPE> ):
			_deapodize_kernel<<< dimGrid, dimBlock >>>( num_elements, imageDevPtr, plan->deapodization_filter_NFFT, fixed_dims_multiplier*plan->domain_size_coils );
			break;
		case sizeof( mr_recon::NFFT_H_plan<UINTd,FLOATd,TYPE> ):
			_deapodize_kernel<<< dimGrid, dimBlock >>>( num_elements, imageDevPtr, plan->deapodization_filter_NFFT_H, fixed_dims_multiplier*plan->domain_size_coils );
			break;
		case sizeof( mr_recon::NFFT_iteration_plan<UINTd,FLOATd,TYPE> ):
			_deapodize_kernel<<< dimGrid, dimBlock >>>( num_elements, imageDevPtr, plan->deapodization_filter_NFFT_it, fixed_dims_multiplier*plan->domain_size_coils );
			break;
		default:
		  printf("Implementation error (deapodize)!!!." );
			break;
		}
	
	}
	else{

		// Cumbersome case: Each image does not match the deapodization filter in size

		// If the last dimensions are fixed then the dimensionality of deapodization filter is reduced
		unsigned int fixed_dims_multiplier = 1, num_elements = 1, num_dims = 0;
		for( unsigned int _d=0; _d<plan->d; _d++ ){
			if( ((unsigned int*)&plan->fixed_dims)[_d] ){
				assert( ((unsigned int*)&plan->matrix_size)[_d] == ((unsigned int*)&plan->matrix_size_os)[_d] );
				fixed_dims_multiplier *= ((unsigned int*)&plan->matrix_size)[_d];
			}
			else{
				num_elements *= ((unsigned int*)&plan->matrix_size_os)[_d];
				num_dims++;
			}
		}

		assert( num_dims<=plan->d );

		dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
		dim3 dimGrid( (unsigned int) ceil((double)num_elements/deviceProp.maxThreadsPerBlock), 1, 1 );

		switch(num_dims){

		case 2:
		  {
		    uint2 filter_size = uintd_to_uint2(plan->matrix_size);
		    uint2 image_size = uintd_to_uint2(plan->matrix_size_os);
		    uint2 corner1 = (image_size-filter_size)>>1;
		    uint2 corner2 = image_size-corner1;

		    switch( sizeof(*plan))
		      {
		      case sizeof( mr_recon::NFFT_plan<UINTd,FLOATd,TYPE> ):
			deapodize_kernel<uint2><<< dimGrid, dimBlock >>> ( image_size, filter_size, corner1, corner2, imageDevPtr, plan->deapodization_filter_NFFT, fixed_dims_multiplier*plan->domain_size_coils );
		      break;
		      case sizeof( mr_recon::NFFT_H_plan<UINTd,FLOATd,TYPE> ):
			deapodize_kernel<uint2><<< dimGrid, dimBlock >>> ( image_size, filter_size, corner1, corner2, imageDevPtr, plan->deapodization_filter_NFFT_H, fixed_dims_multiplier*plan->domain_size_coils );
		      break;
		      case sizeof( mr_recon::NFFT_iteration_plan<UINTd,FLOATd,TYPE> ):
			deapodize_kernel<uint2><<< dimGrid, dimBlock >>> ( image_size, filter_size, corner1, corner2, imageDevPtr, plan->deapodization_filter_NFFT_it, fixed_dims_multiplier*plan->domain_size_coils );
		      break;
		      default:
			printf("Implementation error (deapodize)!!!." );
			break;
		      }
		    break;
		  }
		  
		case 3:
		  {
		    uint3 filter_size = uintd_to_uint3(plan->matrix_size);
		    uint3 image_size = uintd_to_uint3(plan->matrix_size_os);
		    uint3 corner1 = (image_size-filter_size)>>1;
		    uint3 corner2 = image_size-corner1;

		    switch( sizeof(*plan))
		      {
		      case sizeof( mr_recon::NFFT_plan<UINTd,FLOATd,TYPE> ):
			deapodize_kernel<uint3><<< dimGrid, dimBlock >>> ( image_size, filter_size, corner1, corner2, imageDevPtr, plan->deapodization_filter_NFFT, fixed_dims_multiplier*plan->domain_size_coils );
		      break;
		      case sizeof( mr_recon::NFFT_H_plan<UINTd,FLOATd,TYPE> ):
			deapodize_kernel<uint3><<< dimGrid, dimBlock >>> ( image_size, filter_size, corner1, corner2, imageDevPtr, plan->deapodization_filter_NFFT_H, fixed_dims_multiplier*plan->domain_size_coils );
		      break;
		      case sizeof( mr_recon::NFFT_iteration_plan<UINTd,FLOATd,TYPE> ):
			deapodize_kernel<uint3><<< dimGrid, dimBlock >>> ( image_size, filter_size, corner1, corner2, imageDevPtr, plan->deapodization_filter_NFFT_it, fixed_dims_multiplier*plan->domain_size_coils );
		      break;
		      default:
			printf("Implementation error (deapodize)!!!." );
			break;
		      }
		    break;
		  }

		case 4:
		  {
		    uint4 filter_size = uintd_to_uint4(plan->matrix_size);
		    uint4 image_size = uintd_to_uint4(plan->matrix_size_os);
		    uint4 corner1 = (image_size-filter_size)>>1;
		    uint4 corner2 = image_size-corner1;

		    switch( sizeof(*plan))
		      {
		      case sizeof( mr_recon::NFFT_plan<UINTd,FLOATd,TYPE> ):
			deapodize_kernel<uint4><<< dimGrid, dimBlock >>> ( image_size, filter_size, corner1, corner2, imageDevPtr, plan->deapodization_filter_NFFT, fixed_dims_multiplier*plan->domain_size_coils );
		      break;
		      case sizeof( mr_recon::NFFT_H_plan<UINTd,FLOATd,TYPE> ):
			deapodize_kernel<uint4><<< dimGrid, dimBlock >>> ( image_size, filter_size, corner1, corner2, imageDevPtr, plan->deapodization_filter_NFFT_H, fixed_dims_multiplier*plan->domain_size_coils );
		      break;
		      case sizeof( mr_recon::NFFT_iteration_plan<UINTd,FLOATd,TYPE> ):
			deapodize_kernel<uint4><<< dimGrid, dimBlock >>> ( image_size, filter_size, corner1, corner2, imageDevPtr, plan->deapodization_filter_NFFT_it, fixed_dims_multiplier*plan->domain_size_coils );
		      break;
		      default:
			printf("Implementation error (deapodize)!!!." );
			break;
		      }
		    break;
		  }
		}
		
		cudaError_t err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected in 'deapodize_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}

	}
}

// Get density compensation weights
template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ float* get_density_compensation_weights( PLAN<UINTd, FLOATd, TYPE> *plan )
{
	return plan->weights_DevPtr;
}

// Get sample trajectories
template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ FLOATd* get_sample_trajectories( PLAN<UINTd, FLOATd, TYPE> *plan )
{
	FLOATd *result = 0x0;
	cudaMalloc((void**)&result, plan->samples_per_projection*plan->projections_per_frame*sizeof(FLOATd));

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'get_sample_trajectories' (malloc?): %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Compute trajectories
	bool success = NFFT_get_samples( plan, result );

	if( !success ){
		cudaFree(result);
		return 0x0;
	}
	else
		return result;
}


/*
	Class implementation of the initialize, convolve, compute, and cleanup methods for the three plans supported.
*/

// Initialize NFFT

template< class UINTd, class FLOATd, char TYPE > __host__ bool
mr_recon::NFFT_plan< UINTd, FLOATd, TYPE >::initialize()
{	
	printf("\nNFFT_plan::initialize() not yet implemented.\n");
	return false;
}

__host__ float
computeBeta( unsigned int matrix_size, unsigned int matrix_size_os, float W )
{	
  // Compute Kaiser-Bessel beta paramter according to the formula provided in 
  // Beatty et. al. IEEE TMI 2005;24(6):799-808.
  double alpha = (double)matrix_size_os/(double)matrix_size;
  double _beta = (CUDART_PI*sqrt((((double)W)*((double)W))/(alpha*alpha)*(alpha-0.5)*(alpha-0.5)-0.8)); 
  //double _beta = CUDART_PI*sqrt((((double)W)*((double)(W))*(alpha-0.5)*(alpha-0.5))/(alpha*alpha)-1.0);
//  printf("\nSetting beta to %f", _beta);
  return (float) _beta;
}


// Initialize NFFT_H

template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_H_plan< UINTd, FLOATd, TYPE >::initialize()
{	
	initialized = false;

	beta = computeBeta( matrix_size.x, matrix_size_os.x, W );

	// Calculate deapodization filter
	if( calculate_deapodization ){

		// If the last dimensions are fixed then the dimensionality of deapodization filter is reduced
		unsigned int num_elements = 1;
		for( unsigned int _d=0; _d<d; _d++ ){
			if( ((unsigned int*)&fixed_dims)[_d]==0 )
				num_elements *= ((unsigned int*)&matrix_size)[_d];
		}

		cudaMalloc((void**)&deapodization_filter_NFFT_H, num_elements*sizeof(cuFloatComplex));
		cuda_calculate_deapodization_filter( matrix_size, matrix_size_os, fixed_dims, W, domain_size_grid, deapodization_filter_NFFT_H );

		//DEBUG
//		printf("\nWriting deapod.raw with %d elements.\n", num_elements );
//		float* tmp = (float*) calloc( 1, prod(matrix_size)*sizeof(cuFloatComplex) );
//		cudaMemcpy( tmp, deapodization_filter_NFFT_H, prod(matrix_size)*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
//		FILE *fout = fopen("deapod.raw", "wb");
//		fwrite( tmp, prod(matrix_size), 2*sizeof(float), fout );
//		fclose(fout);
	}

	initialized = true;
	return true;
}


// Initialize NFFT iteration

template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_iteration_plan< UINTd, FLOATd, TYPE >::initialize()
{	
	initialized = false;

	beta = computeBeta( matrix_size.x, matrix_size_os.x, W );

	// Calculate deapodization filter (if the last dimensions are fixed then the dimensionality of deapodization filter is reduced).
	unsigned int num_elements = 1;
	for( unsigned int _d=0; _d<d; _d++ ){
		if( ((unsigned int*)&fixed_dims)[_d]==0 )
			num_elements *= ((unsigned int*)&matrix_size)[_d];
	}

	cudaMalloc((void**)&deapodization_filter_NFFT_it, num_elements*sizeof(cuFloatComplex));
	cuda_calculate_deapodization_filter( matrix_size, matrix_size_os, fixed_dims, W, domain_size_grid, deapodization_filter_NFFT_it );
/*
	//DEBUG
	printf("\nWriting deapod.raw with %d elements.\n", num_elements );
	float* tmp = (float*) calloc( 1, num_elements*sizeof(cuFloatComplex) );
	cudaMemcpy( tmp, deapodization_filter_NFFT_it, num_elements*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
	for(int i=0; i<num_elements; i++)
	  tmp[i]=sqrtf(tmp[2*i]*tmp[2*i]+tmp[2*i+1]*tmp[2*i+1]);
	FILE *fout = fopen("deapod.raw", "wb");
	fwrite( tmp, num_elements, sizeof(float), fout );
	fclose(fout);
	free(tmp);
*/

	initialized = true;
	return true;
}

// Compute NFFT

template <class UINTd, class FLOATd, char TYPE> __host__ bool
mr_recon::NFFT_plan< UINTd, FLOATd, TYPE >::compute( cuFloatComplex *samplesDevPtr, cuFloatComplex *imageDevPtr, bool oversampled_image, bool densityCompensate )
{
	printf("\nNFFT_plan::compute() not yet implemented (only in the iteration plan)\n");
	return false;
}


// Compute NFFT_H
template <class UINTd, class FLOATd, char TYPE> __host__ bool
mr_recon::NFFT_H_plan< UINTd, FLOATd, TYPE >::compute( cuFloatComplex *samplesDevPtr, cuFloatComplex *_imageDevPtr, bool oversampled_image, bool densityCompensate )
{
	// Sanity check
	for( unsigned int i=0; i<d; i++ ){
		if( ((unsigned int*)&fixed_dims)[i] )
			if( ((unsigned int*) &domain_size_grid)[i] > 1 ){
				printf("\nConvolution error: Cannot have fixed dim with corresponding domain dimension > 1. Ignoring.");
				return false;
			}
	}

	cuFloatComplex *imageDevPtr;

	if( !oversampled_image ){
		cudaMalloc((void**)&imageDevPtr, domain_size_coils*prod(matrix_size_os)*sizeof(cuFloatComplex));
	}
	else
		imageDevPtr = _imageDevPtr;

	// Density compensation
	if( densityCompensate )
		if( weights_DevPtr )
			density_compensate( this, samplesDevPtr, samplesDevPtr );
		else
			printf("\nWARNING: densityCompensation requested but no weights have been uploaded or computed. Ignoring.");

#ifdef VISUAL_TRACE
			printf("\nBefore convolution (H)"); fflush(stdout);
#endif
	// Convolution
	bool success = NFFT_convolve( this, samplesDevPtr, imageDevPtr );

#ifdef VISUAL_TRACE
			printf("\nAfter convolution (H)"); fflush(stdout);
#endif

	// FFT (if we are currently calculating the deapodization filter we want no FFT scaling).
	cuFloatComplex *tmp;
	cudaMalloc( (void**)&tmp, domain_size_coils*prod(matrix_size_os)*sizeof(cuFloatComplex) );

#ifdef VISUAL_TRACE
	printf("\nBefore FFT"); fflush(stdout);
#endif
/*
	//DEBUG
	float* _tmp = (float*) calloc( 1, domain_size_coils*prod(matrix_size_os)*sizeof(cuFloatComplex) );
	cudaMemcpy( _tmp, imageDevPtr, domain_size_coils*prod(matrix_size_os)*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
	FILE *fout = fopen("GPU_kspace.raw", "wb");
	fwrite( _tmp, domain_size_coils*prod(matrix_size_os), sizeof(cuFloatComplex), fout );
	fclose(fout);
*/
	fft_shift( imageDevPtr, tmp, matrix_size_os, domain_size_coils);
	if( success )
		success = K2I_ALL( tmp, matrix_size_os, domain_size_coils, do_deapodization, false );	
	fft_shift( tmp, imageDevPtr, matrix_size_os, domain_size_coils);
	cudaFree(tmp);

#ifdef VISUAL_TRACE
	printf("\nAfter FFT"); fflush(stdout);
#endif

	// Deapodization (and cropping for non-oversampled output).
	if( oversampled_image ){

		// Deapodization (only when computing the filter itself does the 'do_deapodization if' fail).
		if( success && do_deapodization )
			deapodize( this, true, imageDevPtr );
	}
	else{

		// Crop
		if( success )
			crop_image( matrix_size, matrix_size_os, _imageDevPtr, imageDevPtr, domain_size_coils );

		// Deapodize
		if( success && do_deapodization )
			deapodize( this, false, _imageDevPtr );
		
		cudaFree( imageDevPtr );
	}

	return success;
}


// Compute NFFT_iteration

template <class UINTd, class FLOATd, char TYPE> __host__ bool
mr_recon::NFFT_iteration_plan< UINTd, FLOATd, TYPE >::compute( cuFloatComplex *samplesDevPtr, cuFloatComplex *imageDevPtr, bool oversampled_image, bool densityCompensate )
{

	// Sanity check
	for( unsigned int i=0; i<d; i++ ){
		if( ((unsigned int*)&fixed_dims)[i] )
			if( ((unsigned int*) &domain_size_grid)[i] > 1 ){
				printf("\nConvolution error: Cannot have fixed dim with corresponding domain dimension > 1. Ignoring.");
				return false;
			}
	}

	bool success = true;

	// Deapodization
	deapodize( this, oversampled_image, imageDevPtr );

	// Setup oversampled image (if not given as input)
	cuFloatComplex *image_os;
	if( oversampled_image ){
		image_os = imageDevPtr;
		cuda_border_fill_image<cuFloatComplex,UINTd>( matrix_size, matrix_size_os, make_float2( 0.0f, 0.0f ), image_os, domain_size_coils );
	}
	else{
		cudaMalloc( (void**) &image_os, domain_size_coils*prod(matrix_size_os)*sizeof(cuFloatComplex) );
		image_copy_with_zero_fill( matrix_size, matrix_size_os, imageDevPtr, image_os, domain_size_coils );
	}

	// FFT
	cuFloatComplex *tmp;
	cudaMalloc( (void**)&tmp, domain_size_coils*prod(matrix_size_os)*sizeof(cuFloatComplex) );

#ifdef VISUAL_TRACE
	printf("\nBefore FFT"); fflush(stdout);
#endif

	fft_shift( image_os, tmp, matrix_size_os, domain_size_coils);
	if( success )
		success = I2K_ALL( tmp, matrix_size_os, domain_size_coils, false, false );
	fft_shift( tmp, image_os, matrix_size_os, domain_size_coils);

#ifdef VISUAL_TRACE
	printf("\nAfter FFT"); fflush(stdout);
#endif

#ifdef VISUAL_TRACE
	printf("\nBefore convolution"); fflush(stdout);
#endif

	// Convolution
	if( success )
		success = NFFT_iteration_convolve_to_samples( this, samplesDevPtr, image_os );

#ifdef VISUAL_TRACE
	printf("\nAfter convolution"); fflush(stdout);
#endif

	// Density compensation
	if( densityCompensate )
		if( weights_DevPtr )
			density_compensate( this, samplesDevPtr, samplesDevPtr );
		else
			printf("\nWARNING: densityCompensation requested but no weights have been uploaded or computed. Ignoring.");

	// Convolution

#ifdef VISUAL_TRACE
	printf("\nBefore convolution (H)"); fflush(stdout);
#endif

	if( success )
		success = NFFT_iteration_convolve_to_image( this, samplesDevPtr, image_os );

#ifdef VISUAL_TRACE
	printf("\nAfter convolution (H)"); fflush(stdout);
#endif
#ifdef VISUAL_TRACE
	printf("\nBefore FFT"); fflush(stdout);
#endif

	// FFT
	fft_shift( image_os, tmp, matrix_size_os, domain_size_coils);
	if( success )
		success = K2I_ALL( tmp, matrix_size_os, domain_size_coils, false, false );
	fft_shift( tmp, image_os, matrix_size_os, domain_size_coils);
	cudaFree(tmp);

#ifdef VISUAL_TRACE
	printf("\nAfter FFT"); fflush(stdout);
#endif

	if( oversampled_image ){
	  
		// Deapodization
		if( success ) 
			deapodize( this, true, imageDevPtr );

		// Zero fill
		if( success )
			cuda_border_fill_image( matrix_size, matrix_size_os, make_cuFloatComplex(0.0f, 0.0f), imageDevPtr, domain_size_coils );
	}
	else{

		// Crop
		if( success )
			crop_image( matrix_size, matrix_size_os, imageDevPtr, image_os, domain_size_coils );

		// Deapodization
		if( success ) 
			deapodize( this, false, imageDevPtr );

		cudaFree( image_os );
	}

	return success;
}


// Free all implementation-allocated cuda memory: NFFT

template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_plan< UINTd, FLOATd, TYPE >::cleanup()
{
	// Free all memory

	/* MSH 2010.01.07 RemGlobals
	if( (void*)::sample_positions_f2_NFFT == (void*)sample_positions_DevPtr )
		::sample_positions_f2_NFFT = 0x0;

        

	if( (void*)::sample_positions_f4_NFFT == (void*)sample_positions_DevPtr )
		::sample_positions_f4_NFFT = 0x0;


//	if( ::domainsMap_NFFT == domainsMap_NFFT ) NO DEVICE_PTR CURRENTLY
//		::domainsMap_NFFT = 0x0;

	if( ::stripsMap_NFFT == stripsMapDevPtr_NFFT )
		::stripsMap_NFFT = 0x0;


	if( (void*)::stripsDir_ui2_NFFT == (void*)stripsDirDevPtr_NFFT )
		::stripsDir_ui2_NFFT = 0x0;

	if( (void*)::stripsDir_ui4_NFFT == (void*)stripsDirDevPtr_NFFT )
		::stripsDir_ui4_NFFT = 0x0;

	if( (void*)::stripOrigins_ui2_NFFT == (void*)stripOriginsDevPtr_NFFT )
		::stripOrigins_ui2_NFFT = 0x0;

	if( (void*)::stripOrigins_ui4_NFFT == (void*)stripOriginsDevPtr_NFFT )
		::stripOrigins_ui4_NFFT = 0x0;

	if( ::stripLengths_NFFT == stripLengthsDevPtr_NFFT )
		::stripLengths_NFFT = 0x0;

	*/

	if( sample_positions_DevPtr )
		cudaFree( sample_positions_DevPtr );
	sample_positions_DevPtr = 0x0;

	if( sample_positions )
		delete sample_positions;
	sample_positions = 0x0;

	if( stripsMapDevPtr_NFFT )
		cudaFree( stripsMapDevPtr_NFFT );
	stripsMapDevPtr_NFFT = 0x0;
	
	if( stripsMap_NFFT )
		delete stripsMap_NFFT;
	stripsMap_NFFT = 0x0;

	if( stripsDirDevPtr_NFFT )
		cudaFree( stripsDirDevPtr_NFFT );
	stripsDirDevPtr_NFFT = 0x0;
		
	if( stripsDir_NFFT )
		delete stripsDir_NFFT;
	stripsDir_NFFT = 0x0;

	if( stripOriginsDevPtr_NFFT )
		cudaFree( stripOriginsDevPtr_NFFT );
	stripOriginsDevPtr_NFFT = 0x0;
		
	if( stripOrigins_NFFT )
		delete stripOrigins_NFFT;
	stripOrigins_NFFT = 0x0;

	if( stripLengthsDevPtr_NFFT )
		cudaFree( stripLengthsDevPtr_NFFT );
	stripLengthsDevPtr_NFFT = 0x0;
		
	if( stripLengths_NFFT )
		delete stripLengths_NFFT;
	stripLengths_NFFT = 0x0;

	// TODO: the deapod. filter is GLOBAL! We should have one per plan instead.
	if( deapodization_filter_NFFT )
		cudaFree(deapodization_filter_NFFT);
	deapodization_filter_NFFT = 0x0;

	return true;
}


// Free all implementation-allocated cuda memory: NFFT_H

template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_H_plan< UINTd, FLOATd, TYPE >::cleanup()
{
	/*
	if( (void*)::sample_positions_f2_NFFT_H == (void*)sample_positions_DevPtr )
		::sample_positions_f2_NFFT_H = 0x0;

	if( (void*)::sample_positions_f4_NFFT_H == (void*)sample_positions_DevPtr )
		::sample_positions_f4_NFFT_H = 0x0;

	if( ::domainsMap_NFFT_H == domainsMapDevPtr_NFFT_H )
		::domainsMap_NFFT_H = 0x0;

	if( ::stripsMap_NFFT_H == stripsMapDevPtr_NFFT_H )
		::stripsMap_NFFT_H = 0x0;

	if( ::strips_NFFT_H == stripsDevPtr_NFFT_H )
		::strips_NFFT_H = 0x0;

	*/
	if( sample_positions_DevPtr )
		cudaFree( sample_positions_DevPtr );
	sample_positions_DevPtr = 0x0;

	if( sample_positions )
		delete sample_positions;
	sample_positions = 0x0;

	if( weights_DevPtr )
		cudaFree( weights_DevPtr);
	weights_DevPtr = 0x0;

	if( domainsMapDevPtr_NFFT_H )
		cudaFree( domainsMapDevPtr_NFFT_H );
	domainsMapDevPtr_NFFT_H = 0x0;

	if( stripsMapDevPtr_NFFT_H )
		cudaFree( stripsMapDevPtr_NFFT_H );
	stripsMapDevPtr_NFFT_H = 0x0;

	if( stripsMap_NFFT_H )
		delete stripsMap_NFFT_H;
	stripsMap_NFFT_H = 0x0;

	if( stripsDevPtr_NFFT_H )
		cudaFree( stripsDevPtr_NFFT_H );
	stripsDevPtr_NFFT_H = 0x0;

	if( strips_NFFT_H )
		delete strips_NFFT_H;
	strips_NFFT_H = 0x0;

	if( domainsMap_NFFT_H )
		delete domainsMap_NFFT_H;
	domainsMap_NFFT_H = 0x0;

	if( clean_deapodization && deapodization_filter_NFFT_H )
		cudaFree(deapodization_filter_NFFT_H);
	deapodization_filter_NFFT_H = 0x0;

	return true;
}


// Free all implementation-allocated cuda memory: NFFT_iteration

template< class UINTd, class FLOATd, char TYPE > bool
mr_recon::NFFT_iteration_plan< UINTd, FLOATd, TYPE >::cleanup()
{
	// Free all memory

	/*
	if( (void*)::sample_positions_f2_it == (void*)sample_positions_DevPtr )
		::sample_positions_f2_it = 0x0;

	if( (void*)::sample_positions_f4_it == (void*)sample_positions_DevPtr )
		::sample_positions_f4_it = 0x0;

	if( ::stripsMap_NFFT_it == stripsMapDevPtr_NFFT )
		::stripsMap_NFFT_it = 0x0;

	if( (void*)::stripsDir_ui2_NFFT_it == (void*)stripsDirDevPtr_NFFT )
		::stripsDir_ui2_NFFT_it = 0x0;

	if( (void*)::stripsDir_ui4_NFFT_it == (void*)stripsDirDevPtr_NFFT )
		::stripsDir_ui4_NFFT_it = 0x0;

	if( (void*)::stripOrigins_ui2_NFFT_it == (void*)stripOriginsDevPtr_NFFT )
		::stripOrigins_ui2_NFFT_it = 0x0;

	if( (void*)::stripOrigins_ui4_NFFT_it == (void*)stripOriginsDevPtr_NFFT )
		::stripOrigins_ui4_NFFT_it = 0x0;

	if( ::stripLengths_NFFT_it == stripLengthsDevPtr_NFFT )
		::stripLengths_NFFT_it = 0x0;
	*/

	if( sample_positions_DevPtr )
		cudaFree( sample_positions_DevPtr );
	sample_positions_DevPtr = 0x0;

	if( sample_positions )
		delete sample_positions;
	sample_positions = 0x0;

	if( weights_DevPtr )
		cudaFree( weights_DevPtr);
	weights_DevPtr = 0x0;

	if( stripsMapDevPtr_NFFT )
		cudaFree( stripsMapDevPtr_NFFT );
	stripsMapDevPtr_NFFT = 0x0;
	
	if( stripsMap_NFFT )
		delete stripsMap_NFFT;
	stripsMap_NFFT = 0x0;

	if( stripsDirDevPtr_NFFT )
		cudaFree( stripsDirDevPtr_NFFT );
	stripsDirDevPtr_NFFT = 0x0;
		
	if( stripsDir_NFFT )
		delete stripsDir_NFFT;
	stripsDir_NFFT = 0x0;

	if( stripOriginsDevPtr_NFFT )
		cudaFree( stripOriginsDevPtr_NFFT );
	stripOriginsDevPtr_NFFT = 0x0;
		
	if( stripOrigins_NFFT )
		delete stripOrigins_NFFT;
	stripOrigins_NFFT = 0x0;

	if( stripLengthsDevPtr_NFFT )
		cudaFree( stripLengthsDevPtr_NFFT );
	stripLengthsDevPtr_NFFT = 0x0;
		
	if( stripLengths_NFFT )
		delete stripLengths_NFFT;
	stripLengths_NFFT = 0x0;

	/*
	if( ::domainsMap_NFFT_H_it == domainsMapDevPtr_NFFT_H )
		::domainsMap_NFFT_H_it = 0x0;

	if( ::stripsMap_NFFT_H_it == stripsMapDevPtr_NFFT_H )
		::stripsMap_NFFT_H_it = 0x0;

	if( ::strips_NFFT_H_it == stripsDevPtr_NFFT_H )
		::strips_NFFT_H_it = 0x0;
	*/

	if( sample_positions_DevPtr )
		cudaFree( sample_positions_DevPtr );
	sample_positions_DevPtr = 0x0;

	if( sample_positions )
		delete sample_positions;
	sample_positions = 0x0;

	if( weights_DevPtr )
		cudaFree( weights_DevPtr);
	weights_DevPtr = 0x0;

	if( domainsMapDevPtr_NFFT_H )
		cudaFree( domainsMapDevPtr_NFFT_H );
	domainsMapDevPtr_NFFT_H = 0x0;

	if( stripsMapDevPtr_NFFT_H )
		cudaFree( stripsMapDevPtr_NFFT_H );
	stripsMapDevPtr_NFFT_H = 0x0;

	if( stripsMap_NFFT_H )
		delete stripsMap_NFFT_H;
	stripsMap_NFFT_H = 0x0;

	if( stripsDevPtr_NFFT_H )
		cudaFree( stripsDevPtr_NFFT_H );
	stripsDevPtr_NFFT_H = 0x0;

	if( strips_NFFT_H )
		delete strips_NFFT_H;
	strips_NFFT_H = 0x0;

	if( domainsMap_NFFT_H )
		delete domainsMap_NFFT_H;
	domainsMap_NFFT_H = 0x0;

	// TODO: the deapod. filter is GLOBAL! We should have one per plan instead.
	if( deapodization_filter_NFFT_it )
		cudaFree(deapodization_filter_NFFT_it);
	deapodization_filter_NFFT_it = 0x0;

	return true;
}


/* 
	Plans' private member function implementations
*/


template <class UINTd, class FLOATd, char TYPE> __host__ bool
mr_recon::NFFT_plan< UINTd, FLOATd, TYPE >::convolve( cuFloatComplex *samplesDevPtr, cuFloatComplex *result )
{
	// TODO: as below but other direction (only iteration NFFT implemented so far)
	return false;
}

template <class UINTd, class FLOATd, char TYPE> __host__ bool
mr_recon::NFFT_H_plan< UINTd, FLOATd, TYPE >::convolve( cuFloatComplex *samplesDevPtr, cuFloatComplex *imageDevPtr )
{
	return NFFT_convolve_to_image<UINTd, FLOATd, TYPE, true>( this, &tex_sample_values_NFFT_H, samplesDevPtr, imageDevPtr );
}

template <class UINTd, class FLOATd, char TYPE> __host__ bool
mr_recon::NFFT_iteration_plan< UINTd, FLOATd, TYPE >::convolve( cuFloatComplex *samplesDevPtr, cuFloatComplex *result )
{
	printf("\nWARNING: \"NFFT_convolve( <NFFT_iteration_plan>, ... )\" is not well defined. Which direction to go?");
	printf("\nIgnoring call (use 'NFFT_iteration_convolve_to_samples' or 'NFFT_iteration_convolve_to_image' instead.\n");

	return false;
}


// Invoke convolution (NFFT)

template< class UINTd, class FLOATd, char TYPE, bool CORE, template< class, class, char > class PLAN > __host__ bool 
NFFT_convolve_to_samples( PLAN<UINTd, FLOATd, TYPE> *plan, texture<cuFloatComplex, 1, cudaReadModeElementType> *tex_image_values, cuFloatComplex *samplesDevPtr, cuFloatComplex *imageDevPtr )
{
	// Get device properties
	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	const unsigned int warp_size = deviceProp.warpSize;

	/*
		Setup grid and threads
	*/

	// Determined from the .cubin
	size_t shared_mem_overhead;

	// Determined from the Cuda Occupancy Calculator (setting the shared memory use to 0)
	size_t max_threads_per_block;

	switch( sizeof(UINTd) )
	{
	case sizeof(uint2):
		shared_mem_overhead = NFFT_SHARED_MEMORY_OVERHEAD_PER_2D_KERNEL_BLOCK;
		if( (plan->domain_size_coils>1) && (plan->domain_size_samples>1) )
			max_threads_per_block = NFFT_THREADS_PER_2D_KERNEL_GEN;
		else if( plan->domain_size_samples == 1 )
			max_threads_per_block = NFFT_THREADS_PER_2D_KERNEL_PAR;
		else
			max_threads_per_block = NFFT_THREADS_PER_2D_KERNEL_SEQ;
		break;
	case sizeof(uint3):
		shared_mem_overhead = NFFT_SHARED_MEMORY_OVERHEAD_PER_3D_KERNEL_BLOCK;
		if( (plan->domain_size_coils>1) && (plan->domain_size_samples>1) )
			max_threads_per_block = NFFT_THREADS_PER_3D_KERNEL_GEN;
		else if( plan->domain_size_samples == 1 )
			if( deviceProp.regsPerBlock == 8192 )
				max_threads_per_block = NFFT_THREADS_PER_3D_KERNEL_PAR_G80;
			else if( deviceProp.regsPerBlock == 16384 )
				max_threads_per_block = NFFT_THREADS_PER_3D_KERNEL_PAR_G200;
			else{
				printf("\nERROR: unsupported device!");
				return false;
			}
		else
			max_threads_per_block = NFFT_THREADS_PER_3D_KERNEL_SEQ;
		break;
	case sizeof(uint4):
		shared_mem_overhead = NFFT_SHARED_MEMORY_OVERHEAD_PER_4D_KERNEL_BLOCK;
		if( (plan->domain_size_coils>1) && (plan->domain_size_samples>1) )
			max_threads_per_block = NFFT_THREADS_PER_4D_KERNEL_GEN;
		else if( plan->domain_size_samples == 1 )
			max_threads_per_block = NFFT_THREADS_PER_4D_KERNEL_PAR;
		else
			max_threads_per_block = NFFT_THREADS_PER_4D_KERNEL_SEQ;
		break;
	default:
		printf("\nDimensionality not supported!\n");
		return false;
	}

	// Calculate how much memory we can use per thread block (a multiple of 512 bytes on the G8x)
	size_t avail_memory_per_block = ((deviceProp.sharedMemPerBlock-shared_mem_overhead)/MIN_SHARED_MEM_PER_BLOCK)*MIN_SHARED_MEM_PER_BLOCK;

	// Specify how much memory each thread consumes
	size_t bytes_per_thread;
	if( plan->domain_size_coils == 1 )
		bytes_per_thread = plan->domain_size_samples * sizeof(cuFloatComplex);
	else if( plan->domain_size_samples==1 )
		bytes_per_thread = plan->domain_size_coils * sizeof(cuFloatComplex);
	else
		bytes_per_thread = plan->domain_size_samples * plan->domain_size_coils * sizeof(cuFloatComplex);

	// Now which will be the limiting factor, shared memory or registers?
	size_t max_domains_mem = avail_memory_per_block / bytes_per_thread;
	size_t max_domains_reg = max_threads_per_block;

	max_threads_per_block = min((int)max_domains_mem, (int)max_domains_reg);
	max_threads_per_block = (max_threads_per_block/warp_size)*warp_size; // Safety. Should already be guaranteed from the lookup.

#ifdef _DEBUG
//	printf("\nNFFT threads per block: %d (%d x warp_size (%d)). Domain: (%d, %d).\n", max_threads_per_block, max_threads_per_block/warp_size, warp_size, plan->domain_size_samples, plan->domain_size_coils ); fflush(stdout);
#endif

	// Grid dimensions	
	unsigned int gridDimX = (unsigned int) ceil((double)plan->domain_count_samples/(double)max_threads_per_block);

	// Some checks to make sure that the setup went well

	if( max_threads_per_block == 0 ){
		printf("\nERROR: Not enough shared memory for the chosen element and domain sizes! NFFT_H convolution failed.");
		return false;
	}
	if( max_threads_per_block > deviceProp.maxThreadsDim[0])
	{
		printf("\nNFFT IMPLEMENTATION ERROR: Thread dimension exceeded.");
		return false;
	}
	if( gridDimX > (unsigned int)deviceProp.maxGridSize[0])
	{
		printf("\nNFFT IMPLEMENTATION ERROR: Thread grid dimension exceeded.");
		return false;
	}

	// Block and Grid dimensions
	dim3 dimBlock( (unsigned int)max_threads_per_block, 1, 1 ); 
	dim3 dimGrid( gridDimX, 1, 1 );

	// Assign configuration values to constant device memory variables 
	unsigned int _warp_size = deviceProp.warpSize;
	const cudaChannelFormatDesc desc_complexf = cudaCreateChannelDesc<cuFloatComplex>();
	
	// Bind textures
	const cudaChannelFormatDesc desc_uint = cudaCreateChannelDesc<unsigned int>();
	const cudaChannelFormatDesc desc_uint2 = cudaCreateChannelDesc<uint2>();
	const cudaChannelFormatDesc desc_uint4 = cudaCreateChannelDesc<uint4>();
	const cudaChannelFormatDesc desc_float2 = cudaCreateChannelDesc<float2>();
	const cudaChannelFormatDesc desc_float4 = cudaCreateChannelDesc<float4>();
/*
	if( TYPE == 0 ){

		if(CORE){
			cudaBindTexture( 0, &tex_stripsMap_NFFT, plan->stripsMap_NFFT, &desc_uint2, plan->domain_count_samples*sizeof(uint2) );
		}
		else
			cudaBindTexture( 0, &tex_stripsMap_NFFT_it, plan->stripsMap_NFFT, &desc_uint2, plan->domain_count_samples*sizeof(uint2) );

		if (CORE){
			cudaBindTexture( 0, &tex_stripLengths_NFFT, plan->stripLengths_NFFT, &desc_uint, plan->number_of_strips_NFFT*sizeof(unsigned int) );
		}
		else
			cudaBindTexture( 0, &tex_stripLengths_NFFT_it, plan->stripLengths_NFFT, &desc_uint, plan->number_of_strips_NFFT*sizeof(unsigned int) );

		//FLOATd **sample_positions_ptr = 0x0;

		switch(sizeof(FLOATd)){
	case sizeof(float2):
		//sample_positions_ptr = (FLOATd**) &sample_positions_f2_it;

		if(CORE){
			cudaBindTexture( 0, &tex_sample_positions_f2_NFFT, (float2*) plan->sample_positions_DevPtr, &desc_float2, plan->number_of_samples*sizeof(float2));
		}
		else
			cudaBindTexture( 0, &tex_sample_positions_f2_it, (float2*) plan->sample_positions_DevPtr, &desc_float2, plan->number_of_samples*sizeof(float2));

		if(CORE){
			cudaBindTexture( 0, &tex_stripsDir_ui2_NFFT, plan->stripsDirDevPtr_NFFT, &desc_uint2, plan->domain_count_samples*sizeof(uint2));
		}
		else 
			cudaBindTexture( 0, &tex_stripsDir_ui2_NFFT_it, plan->stripsDirDevPtr_NFFT, &desc_uint2, plan->domain_count_samples*sizeof(uint2));

		if(CORE){
			cudaBindTexture( 0, &tex_stripOrigins_ui2_NFFT, plan->stripOriginsDevPtr_NFFT, &desc_uint2, plan->number_of_strips_NFFT*sizeof(uint2));
		}
		else
			cudaBindTexture( 0, &tex_stripOrigins_ui2_NFFT_it, plan->stripOriginsDevPtr_NFFT, &desc_uint2, plan->number_of_strips_NFFT*sizeof(uint2));
		break;

	case sizeof(float3):
	case sizeof(float4):
		//sample_positions_ptr = (FLOATd**) &sample_positions_f4_it;

		if(CORE){
			cudaBindTexture( 0, &tex_sample_positions_f4_NFFT, (float4*) plan->sample_positions_DevPtr, &desc_float4, plan->number_of_samples*sizeof(float4));
		}
		else
			cudaBindTexture( 0, &tex_sample_positions_f4_it, (float4*) plan->sample_positions_DevPtr, &desc_float4, plan->number_of_samples*sizeof(float4));

		if(CORE){
			cudaBindTexture( 0, &tex_stripsDir_ui4_NFFT, plan->stripsDirDevPtr_NFFT, &desc_uint4, plan->domain_count_samples*sizeof(uint4));
		}
		else
			cudaBindTexture( 0, &tex_stripsDir_ui4_NFFT_it, plan->stripsDirDevPtr_NFFT, &desc_uint4, plan->domain_count_samples*sizeof(uint4)); 

		if(CORE){
			cudaBindTexture( 0, &tex_stripOrigins_ui4_NFFT, plan->stripOriginsDevPtr_NFFT, &desc_uint4, plan->number_of_strips_NFFT*sizeof(uint4));
		}
		else
			cudaBindTexture( 0, &tex_stripOrigins_ui4_NFFT_it, plan->stripOriginsDevPtr_NFFT, &desc_uint4, plan->number_of_strips_NFFT*sizeof(uint4));
		break;
		}		
	}
*/
	// Bind sample values texture to input device ptr
	cudaBindTexture( 0, tex_image_values, imageDevPtr, &desc_complexf, plan->domain_size_coils*prod(plan->matrix_size_os)*sizeof(cuFloatComplex));

	/*
		Invoke kernel
	*/

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda NFFT error before convolution: %s\n", cudaGetErrorString(err) ); fflush(stdout);
		return false;
	}

	if( (plan->domain_size_coils>1) && (plan->domain_size_samples>1) )
	{
		// Generic kernel
		printf("\nError: disabled convolution.<\n");
//	  NFFT_convolve_kernel_generic<UINTd,FLOATd,0,CORE><<< dimGrid, dimBlock, prod(dimBlock)*bytes_per_thread >>>
//			( _half_alpha_times_W_SQR, _one_over_alpha_squared, _one_over_W, _beta, plan->domain_size_samples, plan->domain_count_samples, plan->domain_size_coils, plan->matrix_size_os, _warp_size, plan->number_of_strips_NFFT, samplesDevPtr );

	  cudaError_t err = cudaGetLastError();
	  if( err != cudaSuccess ){
		  printf("\nCuda error detected in 'NFFT_convolve_kernel_generic': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		  exit(1);
	  }
	}
	else if( plan->domain_size_samples == 1 )
	{
		// Parallel kernel

		if( TYPE == 1 || TYPE == 2 )
			_nfft_NFFT_set_constants( plan );

		unsigned int double_warp_size_power=0;
		unsigned int __tmp = _warp_size<<1;
		while(__tmp!=1){
			__tmp>>=1;
			double_warp_size_power++;
		}

		cudaError_t err;

		err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected before 'NFFT_convolve_kernel_parallel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}

		if( TYPE == 0 )
		        NFFT_convolve_kernel_parallel_generic<UINTd,FLOATd,TYPE,CORE><<< dimGrid, dimBlock, prod(dimBlock)*bytes_per_thread >>>
				( plan->alpha, plan->W, beta, plan->domain_count_samples, plan->domain_size_coils, plan->matrix_size, plan->matrix_size_os, plan->matrix_size_wrap, double_warp_size_power, raw_pointer_cast(&(*plan->traj_positions)[0]), samplesDevPtr, (get_last_dim(plan->fixed_dims)==1) );
		else
			NFFT_convolve_kernel_parallel<UINTd,FLOATd,TYPE,CORE><<< dimGrid, dimBlock, prod(dimBlock)*bytes_per_thread >>>
				( plan->alpha, plan->W, beta, plan->domain_count_samples, plan->domain_size_coils, plan->matrix_size, plan->matrix_size_os, double_warp_size_power, 1.0f/(float)plan->projections_per_frame, plan->number_of_strips_NFFT, plan->stripsMapDevPtr_NFFT, plan->stripsDirDevPtr_NFFT, plan->stripOriginsDevPtr_NFFT, plan->stripLengthsDevPtr_NFFT, samplesDevPtr, (get_last_dim(plan->fixed_dims)==1) );

		err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected in 'NFFT_convolve_kernel_parallel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}
	}
	else
	{
		// Sequential kernel
		printf("\nError: disabled convolution.<\n");
//	  NFFT_convolve_kernel_sequential<UINTd,FLOATd,0,CORE><<< dimGrid, dimBlock, prod(dimBlock)*bytes_per_thread >>>
//			( _half_alpha_times_W_SQR, _one_over_alpha_squared, _one_over_W, _beta, plan->domain_size_samples, plan->domain_count_samples, plan->matrix_size_os, _warp_size, plan->number_of_strips_NFFT, samplesDevPtr );

	  cudaError_t err = cudaGetLastError();
	  if( err != cudaSuccess ){
		  printf("\nCuda error detected in 'NFFT_convolve_kernel_sequential': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		  exit(1);
	  }
	}


#ifdef _DEBUG
	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda convolution error: %s\n", cudaGetErrorString(err) ); fflush(stdout);
		return false;
	}
#endif

	// Unbind samples texture
	cudaUnbindTexture( tex_image_values );
/*
	if( TYPE == 0 ){
		if(CORE){
			cudaUnbindTexture(tex_stripsMap_NFFT); }
		else{
			cudaUnbindTexture(tex_stripsMap_NFFT_it); }
		if(CORE){
			cudaUnbindTexture(tex_stripLengths_NFFT); }
		else{
			cudaUnbindTexture(tex_stripLengths_NFFT_it);}

		switch(sizeof(FLOATd)){
	case sizeof(float2):
		if(CORE){
			cudaUnbindTexture(tex_sample_positions_f2_NFFT); }
		else{
			cudaUnbindTexture(tex_sample_positions_f2_it); }
		if(CORE){
			cudaUnbindTexture(tex_stripsDir_ui2_NFFT); }
		else{
			cudaUnbindTexture(tex_stripsDir_ui2_NFFT_it); }
		if(CORE){
			cudaUnbindTexture(tex_stripOrigins_ui2_NFFT); }
		else{
			cudaUnbindTexture(tex_stripOrigins_ui2_NFFT_it); }
		break;
	case sizeof(float3):
	case sizeof(float4):
		if(CORE){
			cudaUnbindTexture(tex_sample_positions_f4_NFFT); }
		else{
			cudaUnbindTexture(tex_sample_positions_f4_it); }
		if(CORE){
			cudaUnbindTexture(tex_stripsDir_ui4_NFFT); }
		else{
			cudaUnbindTexture(tex_stripsDir_ui4_NFFT_it); }
		if(CORE){
			cudaUnbindTexture(tex_stripOrigins_ui4_NFFT); }
		else{
			cudaUnbindTexture(tex_stripOrigins_ui4_NFFT_it); }
		break;
		}		
	}
*/	
	return true;
}

template <class UINTd, class FLOATd, char TYPE> __global__ void
NFFT_get_samples_kernel( unsigned int number_of_samples, UINTd matrix_size, FLOATd *result )
{
	// Global thread number	
	const unsigned int globalThreadId = (blockIdx.x*blockDim.x+threadIdx.x);

	// We might have a couple of warps left over due to our selection of grid/block dimensions.
	// This check will only work when there is NO __SYNCTHREADS() in the code!
	if( globalThreadId >= number_of_samples )
		return;

	const UINTd bias = matrix_size>>1;

	// Sample position to convolve onto
	FLOATd sample_position;
	switch(TYPE){
		case 1:
		case 2:
			sample_position = cuda_compute_radial_sample_position<UINTd, FLOATd, TYPE>( globalThreadId, 1.0f, bias, bias );
			break;
		default:
			return;
	}

	result[globalThreadId] = sample_position;
}

template< class UINTd, class FLOATd, char TYPE, template< class, class, char > class PLAN > 
__host__ bool NFFT_get_samples( PLAN<UINTd, FLOATd, TYPE> *plan, FLOATd* result )
{
	if( TYPE == 0 ){

		// Easy (already computed) - but not yet supported

		return false;
	}

	// Block and Grid dimensions
	dim3 dimBlock( 256 ); 
	dim3 dimGrid( (unsigned int) ceil((double)plan->number_of_samples/(double)dimBlock.x) );

	/*
		Invoke kernel
	*/

	if( TYPE == 1 || TYPE == 2 )
		_nfft_NFFT_set_constants( plan );

	NFFT_get_samples_kernel<UINTd,FLOATd,TYPE><<< dimGrid, dimBlock >>>
		( plan->number_of_samples, plan->matrix_size, result );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'NFFT_get_samples_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
	
	return true;
}

dim3 computeDimGrid( uint2 dimBlock, uint2 domain_count, uint2 domain_size )
{
	if( (domain_count.x*domain_size.x)%dimBlock.x || (domain_count.y*domain_size.y)%dimBlock.y ){
		printf("\ndFATAL ERROR: dimBlock does not match domain size! Quitting.\n");
		exit(1);
	}

	return dim3( (domain_count.x*domain_size.x)/dimBlock.x, (domain_count.y*domain_size.y)/dimBlock.y, 1 );
}

dim3 computeDimGrid( uint3 dimBlock, uint3 domain_count, uint3 domain_size )
{
	if( (domain_count.x*domain_size.x)%dimBlock.x || (domain_count.y*domain_count.z*domain_size.y*domain_size.z)%(dimBlock.y*dimBlock.z) ){
		printf("\ndFATAL ERROR: dimBlock does not match domain size! Quitting.\n");
		exit(1);
	}

	return dim3( (domain_count.x*domain_size.x)/dimBlock.x, (domain_count.y*domain_count.z*domain_size.y*domain_size.z)/(dimBlock.y*dimBlock.z), 1 );
}

dim3 computeDimGrid( uint4 dimBlock, uint4 domain_count, uint4 domain_size )
{
	if( (domain_count.x*domain_size.x)%dimBlock.x || (domain_count.y*domain_count.z*domain_count.w*domain_size.y*domain_size.z*domain_size.w)%(dimBlock.y*dimBlock.z) ){
		printf("\ndFATAL ERROR: dimBlock does not match domain size! Quitting.\n");
		exit(1);
	}
	return dim3( (domain_count.x*domain_size.x)/dimBlock.x, (domain_count.y*domain_count.z*domain_count.w*domain_size.y*domain_size.z*domain_size.w)/(dimBlock.y*dimBlock.z), 1 );
}

// Invoke convolution (NFFT_H)

template< class UINTd, class FLOATd, char TYPE, bool CORE, template< class, class, char > class PLAN > __host__ bool 
NFFT_convolve_to_image( PLAN<UINTd, FLOATd, TYPE> *plan, texture<cuFloatComplex, 1, cudaReadModeElementType> *sample_values, cuFloatComplex *samplesDevPtr, cuFloatComplex *imageDevPtr, bool accumulate )
{
	// Get device properties
	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	const unsigned int warp_size = deviceProp.warpSize;

	/*
		Setup grid and threads
	*/

	// Determined from the .cubin
	size_t shared_mem_overhead;

	// Determined from the Cuda Occupancy Calculator (setting the shared memory use to 0)
	size_t max_threads_per_block;

	switch( sizeof(UINTd) )
	{
	case sizeof(uint2):
		shared_mem_overhead = NFFT_H_SHARED_MEMORY_OVERHEAD_PER_2D_KERNEL_BLOCK;
		if( (plan->domain_size_coils>1) && (prod(plan->domain_size_grid)>1) ){
			if( plan->domain_size_grid.x == prod(plan->domain_size_grid) )
				max_threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_SEMI_GEN;
			else
				max_threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_GEN;
		}
		else if( prod(plan->domain_size_grid) == 1 )
			max_threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_PAR;
		else if( plan->domain_size_coils == 1 )
			max_threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_SEQ;
		break;
	case sizeof(uint3):
		shared_mem_overhead = NFFT_H_SHARED_MEMORY_OVERHEAD_PER_3D_KERNEL_BLOCK;
		if( (plan->domain_size_coils>1) && (prod(plan->domain_size_grid)>1) ){
			if( plan->domain_size_grid.x == prod(plan->domain_size_grid) )
				max_threads_per_block = NFFT_H_THREADS_PER_3D_KERNEL_SEMI_GEN;
			else
				max_threads_per_block = NFFT_H_THREADS_PER_3D_KERNEL_GEN;
		}
		else if( prod(plan->domain_size_grid) == 1 )
			if( deviceProp.regsPerBlock == 8192 )
				max_threads_per_block = NFFT_H_THREADS_PER_3D_KERNEL_PAR_G80;
			else if( deviceProp.regsPerBlock == 16384 )
				max_threads_per_block = NFFT_H_THREADS_PER_3D_KERNEL_PAR_G200;
			else{
				printf("\nERROR: unsupported device!");
				return false;
			}
		else if( plan->domain_size_coils == 1 )
			max_threads_per_block = NFFT_H_THREADS_PER_3D_KERNEL_SEQ;
		break;
	case sizeof(uint4):
		shared_mem_overhead = NFFT_H_SHARED_MEMORY_OVERHEAD_PER_4D_KERNEL_BLOCK;
		if( (plan->domain_size_coils>1) && (prod(plan->domain_size_grid)>1) ){
			if( plan->domain_size_grid.x == prod(plan->domain_size_grid) )
				max_threads_per_block = NFFT_H_THREADS_PER_4D_KERNEL_SEMI_GEN;
			else
				max_threads_per_block = NFFT_H_THREADS_PER_4D_KERNEL_GEN;
		}
		else if( prod(plan->domain_size_grid) == 1 )
			max_threads_per_block = NFFT_H_THREADS_PER_4D_KERNEL_PAR;
		else if( plan->domain_size_coils == 1 )
			max_threads_per_block = NFFT_H_THREADS_PER_4D_KERNEL_SEQ;
		break;
	default:
		printf("\nDimensionality not supported!\n");
		return false;
	}

	// Calculate how much memory we can use per thread block (a multiple of 512 bytes on the G8x)
	size_t avail_memory_per_block = ((deviceProp.sharedMemPerBlock-shared_mem_overhead)/MIN_SHARED_MEM_PER_BLOCK)*MIN_SHARED_MEM_PER_BLOCK;

	// Specify how much memory each thread (corresponding to a grid cell) consumes
	size_t bytes_per_thread = plan->domain_size_coils * sizeof(cuFloatComplex);

	// Now which will be the limiting factor, shared memory or registers?
	size_t max_domains_mem = avail_memory_per_block / bytes_per_thread;
	size_t max_domains_reg = max_threads_per_block;
	
	max_threads_per_block = min( (int)max_domains_mem, (int)max_domains_reg );
	max_threads_per_block = (max_threads_per_block/warp_size)*warp_size;

#ifdef _DEBUG
//	printf("\nNFFT_H threads per block: %d (%d x warp_size (%d)). Domain(%d, %d).\n", max_threads_per_block, max_threads_per_block/warp_size, warp_size, prod(plan->domain_size_grid), plan->domain_size_coils ); fflush(stdout);
//	printf("\nmax_domains_mem/max_domains_reg: %d/%d", max_domains_mem, max_domains_reg );
#endif

	// Some checks to make sure that the setup went well

	if( max_threads_per_block == 0 ){
		printf("\nERROR: Not enough shared memory for the chosen element and domain sizes! NFFT_H convolution failed.");
		return false;
	}

	// Assign configuration values to constant device memory variables 
	unsigned int _warp_size = deviceProp.warpSize;
	const cudaChannelFormatDesc desc_complexf = cudaCreateChannelDesc<cuFloatComplex>();

	// Bind textures
  
	const cudaChannelFormatDesc desc_uint = cudaCreateChannelDesc<unsigned int>();
	const cudaChannelFormatDesc desc_uint2 = cudaCreateChannelDesc< uint2 >();
	const cudaChannelFormatDesc desc_float2 = cudaCreateChannelDesc< float2 >();
	const cudaChannelFormatDesc desc_float4 = cudaCreateChannelDesc< float4 >();

	unsigned int number_of_domains_NFFT_H = prod(plan->domain_count_grid);
/*
	if( TYPE == 0 ){
		if(CORE){
			cudaBindTexture( 0, &tex_domainsMap_NFFT_H, plan->domainsMap_NFFT_H, &desc_uint, number_of_domains_NFFT_H*sizeof(unsigned int) );
		}
		else
			cudaBindTexture( 0, &tex_domainsMap_NFFT_H_it, plan->domainsMap_NFFT_H, &desc_uint, number_of_domains_NFFT_H*sizeof(unsigned int) );
		if(CORE){
			cudaBindTexture( 0, &tex_stripsMap_NFFT_H, plan->stripsMap_NFFT_H, &desc_uint2, number_of_domains_NFFT_H*sizeof(uint2) );
		}
		else 
			cudaBindTexture( 0, &tex_stripsMap_NFFT_H_it, plan->stripsMap_NFFT_H, &desc_uint2, number_of_domains_NFFT_H*sizeof(uint2) );
		if(CORE){
			cudaBindTexture( 0, &tex_strips_NFFT_H, plan->strips_NFFT_H, &desc_uint2, plan->number_of_strips_NFFT_H*sizeof(uint2) );
		}
		else
			cudaBindTexture( 0, &tex_strips_NFFT_H_it, plan->strips_NFFT_H, &desc_uint2, plan->number_of_strips_NFFT_H*sizeof(uint2) );

		//FLOATd **sample_positions_ptr = 0x0;

		switch(sizeof(FLOATd)){

	case sizeof(float2):
		//sample_positions_ptr = (CORE) ? (FLOATd**) &sample_positions_f2_NFFT_H : (FLOATd**) &sample_positions_f2_it;
		if(CORE){
			cudaBindTexture( 0, &tex_sample_positions_f2_NFFT_H, (float2*) plan->sample_positions_DevPtr, &desc_float2, plan->number_of_samples*sizeof(FLOATd));
		}
		else
			cudaBindTexture( 0, &tex_sample_positions_f2_it, (float2*) plan->sample_positions_DevPtr, &desc_float2, plan->number_of_samples*sizeof(float2));
		break;

	case sizeof(float3):
	case sizeof(float4):
		//sample_positions_ptr = (CORE) ? (FLOATd**) &sample_positions_f4_NFFT_H : (FLOATd**) &sample_positions_f4_it;
		if(CORE){
			cudaBindTexture( 0, &tex_sample_positions_f4_NFFT_H, (float4*) plan->sample_positions_DevPtr, &desc_float4, plan->number_of_samples*sizeof(float4));
		}
		else
			cudaBindTexture( 0, &tex_sample_positions_f4_it, (float4*) plan->sample_positions_DevPtr, &desc_float4, plan->number_of_samples*sizeof(float4));
		break;
		}
	}
*/
	// Bind samples texture
	cudaBindTexture( 0, sample_values, samplesDevPtr, &desc_complexf, plan->domain_size_coils*plan->number_of_samples*sizeof(cuFloatComplex));

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error before convolution: %s\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	/*
		Invoke kernel
	*/

	if( (plan->domain_size_coils>1) && (prod(plan->domain_size_grid)>1) )
	{
		if( plan->domain_size_grid.x == prod(plan->domain_size_grid) ){

			// Semigeneric kernel

			if( TYPE == 1 || TYPE == 2 )
				_nfft_NFFT_set_constants( plan );

			dim3 dimBlock(128,1,1);
			dim3 dimGrid(1,1,1); // TODO

			// Allocate temporary memory buffer that includes wrap zone
			cuFloatComplex *_tmp;
			cudaMalloc((void**)&_tmp, plan->domain_size_coils*prod(plan->matrix_size_os+plan->matrix_size_wrap)*sizeof(cuFloatComplex));

			cudaError_t err = cudaGetLastError();
			if( err != cudaSuccess ){
				printf("\nCuda error detected before 'NFFT_H_convolve_kernel_semigeneric' (malloc?): %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
				exit(1);
			}

			unsigned int domain_size_grid_power=0, __tmp = plan->domain_size_grid.x;
			while(__tmp!=1){
				__tmp>>=1;
				domain_size_grid_power++;
			}

			unsigned int double_warp_size_power=0;
			__tmp = _warp_size<<1;
			while(__tmp!=1){
				__tmp>>=1;
				double_warp_size_power++;
			}

		printf("\nError: disabled convolution.<\n");
//			NFFT_H_convolve_kernel_semigeneric<UINTd,FLOATd,TYPE,CORE><<< dimGrid, dimBlock, prod(dimBlock)*bytes_per_thread >>>
//				( plan->domain_count_grid, plan->matrix_size_os, plan->matrix_size_wrap, plan->matrix_size>>1, (plan->matrix_size_os+plan->matrix_size_wrap)>>1, plan->domain_size_grid.x, domain_size_grid_power, double_warp_size_power, 1.0f/(float)plan->projections_per_frame, plan->alpha, plan->W, beta, plan->domain_size_coils, plan->number_of_samples, plan->number_of_strips_NFFT_H, plan->stripsMapDevPtr_NFFT_H, plan->stripsDevPtr_NFFT_H, _tmp );

			err = cudaGetLastError();
			if( err != cudaSuccess ){
				printf("\nCuda error detected in 'NFFT_H_convolve_kernel_semigeneric': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
				exit(1);
			}

			image_wrap( plan->matrix_size_os, plan->matrix_size_wrap, plan->domain_size_coils, accumulate, _tmp, imageDevPtr );

			cudaFree(_tmp);
		}
		else{

			// Generic kernel
		printf("\nError: disabled convolution.<\n");
//		  NFFT_H_convolve_kernel_generic<UINTd,FLOATd,0,CORE><<< dimGrid, dimBlock, prod(dimBlock)*bytes_per_thread >>>
//				( _half_alpha_times_W_SQR, _one_over_alpha_squared, _one_over_W, _beta, plan->domain_size_grid, plan->domain_count_grid, plan->domain_size_coils, plan->matrix_size_os, plan->matrix_size_os>>1, _warp_size, plan->number_of_strips_NFFT_H, plan->number_of_samples, imageDevPtr, accumulate );

		  cudaError_t err = cudaGetLastError();
		  if( err != cudaSuccess ){
			  printf("\nCuda error detected in 'NFFT_H_convolve_kernel_generic': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			  exit(1);
		  }
}
	}
	else if( prod(plan->domain_size_grid) == 1 )
	{

		if( TYPE == 1 || TYPE == 2 )
			_nfft_NFFT_set_constants( plan );
		
		// Parallel kernel
		
		// Allocate and clear temporary memory buffer that includes wrap zone
		cuFloatComplex *_tmp;
		cudaMalloc((void**)&_tmp, plan->domain_size_coils*prod(plan->matrix_size_os+plan->matrix_size_wrap)*sizeof(cuFloatComplex));
//		clear_image( plan->domain_size_coils*prod(plan->matrix_size_os+plan->matrix_size_wrap), make_cuFloatComplex(0.0f, 0.0f), _tmp );

		unsigned int double_warp_size_power=0, __tmp = _warp_size<<1;
		while(__tmp!=1){
			__tmp>>=1;
			double_warp_size_power++;
		}

		cudaError_t err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected before 'NFFT_H_convolve_kernel_parallel' (malloc?): %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}
/*
		unsigned int num_threads_per_block;
		switch(plan->d){
		case 2:
			num_threads_per_block = NFFT_H_THREADS_PER_2D_KERNEL_PAR;
			break;
		case 3:
			num_threads_per_block = NFFT_H_THREADS_PER_3D_KERNEL_PAR;
			break;
		case 4:
			num_threads_per_block = NFFT_H_THREADS_PER_4D_KERNEL_PAR;
			break;
		}
*/
#ifndef SELF_SCHEDULING
		dim3 blockDim(256);
		if( plan->domain_size_coils == 8 )
			blockDim.x = 192;
		dim3 gridDim((unsigned int) ceil((double)prod(plan->domain_count_grid)/(double)blockDim.x) );

		if( TYPE == 0){
			NFFT_H_convolve_kernel_parallel_generic<UINTd,FLOATd,TYPE,CORE><<< gridDim, blockDim, prod(blockDim)*bytes_per_thread >>>
				( plan->domain_count_grid, plan->matrix_size_os, plan->matrix_size_wrap, plan->matrix_size>>1, (plan->matrix_size_os+plan->matrix_size_wrap)>>1, plan->alpha, plan->W, beta, plan->domain_size_coils, double_warp_size_power, plan->number_of_samples, plan->domainsMapDevPtr_NFFT_H, raw_pointer_cast(&(*plan->traj_positions)[0]), raw_pointer_cast(&(*plan->tuples_last)[0]), raw_pointer_cast(&(*plan->bucket_begin)[0]), raw_pointer_cast(&(*plan->bucket_end)[0]), _tmp, (get_last_dim(plan->fixed_dims)==1) );
		}
		else{
			NFFT_H_convolve_kernel_parallel_radial<UINTd,FLOATd,TYPE,CORE><<< gridDim, blockDim, prod(blockDim)*bytes_per_thread >>>
				( plan->domain_count_grid, plan->matrix_size_os, plan->matrix_size_wrap, plan->matrix_size>>1, (plan->matrix_size_os+plan->matrix_size_wrap)>>1, plan->alpha, plan->W, beta, plan->domain_size_coils, double_warp_size_power, 1.0f/(float)plan->projections_per_frame, plan->number_of_samples, plan->number_of_strips_NFFT_H, plan->stripsMapDevPtr_NFFT_H, plan->stripsDevPtr_NFFT_H, plan->domainsMapDevPtr_NFFT_H, _tmp, (get_last_dim(plan->fixed_dims)==1) );
		}

/*
unsigned int num_elements = plan->domain_size_coils*prod(plan->matrix_size_os+plan->matrix_size_wrap);
float* write = (float*) calloc( 1, num_elements*sizeof(cuFloatComplex) );
cudaMemcpy( write, _tmp, num_elements*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost );
FILE *fout = fopen("conv.raw", "wb");
fwrite( write, num_elements, sizeof(cuFloatComplex), fout );
fclose(fout);
*/

#else

		dim3 blockDim(WARP_SIZE, BLOCKDIM_Y);
		dim3 gridDim(30); // Number of multiprocessors (GTX280)

		unsigned int _globalPoolThreadCount = prod(plan->domain_count_grid);
		cudaMemcpyToSymbol( globalPoolThreadCount, &_globalPoolThreadCount, sizeof(unsigned int) );

		unsigned int _globalPoolNextThread = 0;
		cudaMemcpyToSymbol( globalPoolNextThread, &_globalPoolNextThread, sizeof(unsigned int) );
 
		NFFT_H_convolve_kernel_parallel<UINTd,FLOATd,TYPE,CORE><<< gridDim, blockDim, prod(blockDim)*bytes_per_thread >>>
			( plan->domain_count_grid, plan->matrix_size_os, plan->matrix_size_wrap, plan->matrix_size>>1, (plan->matrix_size_os+plan->matrix_size_wrap)>>1, plan->alpha, plan->W, beta, plan->domain_size_coils, double_warp_size_power, 1.0f/(float)plan->projections_per_frame, plan->number_of_samples, plan->number_of_strips_NFFT_H, plan->stripsMapDevPtr_NFFT_H, plan->stripsDevPtr_NFFT_H, plan->domainsMapDevPtr_NFFT_H, _tmp, (get_last_dim(plan->fixed_dims)==1) );

#endif
//		NFFT_H_convolve_kernel_parallel<UINTd,FLOATd,TYPE,CORE><<< gridDim, blockDim, prod(blockDim)*bytes_per_thread >>>
//			( plan->domain_count_grid, plan->matrix_size_os, plan->matrix_size_wrap, plan->matrix_size>>1, (plan->matrix_size_os+plan->matrix_size_wrap)>>1, plan->alpha, plan->W, beta, plan->domain_size_coils, double_warp_size_power, 1.0f/(float)plan->projections_per_frame, plan->number_of_samples, plan->number_of_strips_NFFT_H, plan->number_of_threads_NFFT_H, plan->domainsMapDevPtr_NFFT_H, plan->stripsMapDevPtr_NFFT_H, plan->stripsDevPtr_NFFT_H, _tmp );

		err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected in 'NFFT_H_convolve_kernel_parallel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}

		image_wrap( plan->matrix_size_os, plan->matrix_size_wrap, plan->domain_size_coils, accumulate, _tmp, imageDevPtr );

		cudaFree(_tmp);
	}
	else if( plan->domain_size_coils == 1 )
	{
	  // Sequential kernel
		printf("\nError: disabled convolution.<\n");
//	  NFFT_H_convolve_kernel_sequential<UINTd,FLOATd,0,CORE><<< dimGrid, dimBlock, prod(dimBlock)*bytes_per_thread >>>
//			( _half_alpha_times_W_SQR, _one_over_alpha_squared, _one_over_W, _beta, plan->domain_size_grid, plan->domain_count_grid, plan->matrix_size_os, plan->matrix_size_os>>1, _warp_size, plan->number_of_strips_NFFT_H, plan->number_of_samples, imageDevPtr, accumulate );

	  cudaError_t err = cudaGetLastError();
	  if( err != cudaSuccess ){
		  printf("\nCuda error detected in 'NFFT_H_convolve_kernel_sequential': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		  exit(1);
	  }
	}
	else{
		printf("\nImplementation error: no business here!\n");
		return false;
	}

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda convolution error: %s\n", cudaGetErrorString(err) ); fflush(stdout);
		return false;
	}

	// Unbind samples texture
	cudaUnbindTexture( sample_values );
/*
	if( TYPE == 0){

		if(CORE){
			cudaUnbindTexture(tex_domainsMap_NFFT_H);
		}
		else
			cudaUnbindTexture(tex_domainsMap_NFFT_H_it);
		if(CORE){
			cudaUnbindTexture(tex_stripsMap_NFFT_H);
		}
		else
			cudaUnbindTexture(tex_stripsMap_NFFT_H_it);
		if(CORE){
			cudaUnbindTexture(tex_strips_NFFT_H);
		}
		else
			cudaUnbindTexture(tex_strips_NFFT_H_it);


		switch(sizeof(FLOATd)){

	case sizeof(float2):
		if(CORE){
			cudaUnbindTexture(tex_sample_positions_f2_NFFT_H);
		}
		else
			cudaUnbindTexture(tex_sample_positions_f2_it);
		break;

	case sizeof(float3):
	case sizeof(float4):
		if(CORE){
			cudaUnbindTexture(tex_sample_positions_f4_NFFT_H);
		}
		else
			cudaUnbindTexture( tex_sample_positions_f4_it);
		break;
		}
	}
*/
	return true;
}


/*
  Function to calculate the deapodization filter
*/ 

template< class UINTd > __host__ bool 
cuda_calculate_deapodization_filter( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, float W, UINTd domain_size_grid, cuFloatComplex *out_filterDevPtr )
{
	
	// If the last dimensions are fixed then the dimensionality of deapodization filter is reduced
	unsigned int num_dims = 0;
	for( unsigned int _d=0; _d<sizeof(matrix_size)/sizeof(unsigned int); _d++ ){
		if( ((unsigned int*)&fixed_dims)[_d]==0 )
			num_dims++;
	}

	// Report "result"
	bool success = true;


	//
	//	Grid fictitious trajectory with a single sample at the origin
	//

	switch( num_dims ){
	  
	case 2:
			{
				uint2 _matrix_size = uintd_to_uint2(matrix_size);
				uint2 _matrix_size_os = uintd_to_uint2(matrix_size_os);
				uint2 _domain_size_grid = uintd_to_uint2(domain_size_grid);

				// Grid oversampled image
				cuFloatComplex *filter_os;
				cudaMalloc((void**)&filter_os, prod(_matrix_size_os)*sizeof(cuFloatComplex));

				compute_deapodization_filter<uint2,float2>( _matrix_size_os, W, 1.0f/W, beta, filter_os );

				// FFT
				if( success )
					success = K2I_ALL( filter_os, _matrix_size_os, 1, false, true );

				// Crop
				if( success )
					crop_image( _matrix_size, _matrix_size_os, out_filterDevPtr, filter_os );

				cudaFree(filter_os);

				break;
			}
		case 3:
			{
				uint3 _matrix_size = uintd_to_uint3(matrix_size);
				uint3 _matrix_size_os = uintd_to_uint3(matrix_size_os);
				uint3 _domain_size_grid = uintd_to_uint3(domain_size_grid);

				// Grid oversampled image
				cuFloatComplex *filter_os;
				cudaMalloc((void**)&filter_os, prod(_matrix_size_os)*sizeof(cuFloatComplex));

				compute_deapodization_filter<uint3,float3>( _matrix_size_os, W, 1.0f/W, beta, filter_os );

				// FFT
				if( success )
					success = K2I_ALL( filter_os, _matrix_size_os, 1, false, true );

				// Crop
				if( success )
					crop_image( _matrix_size, _matrix_size_os, out_filterDevPtr, filter_os );
				
				cudaFree(filter_os);

				break;
			}
		case 4:
			{
				uint4 _matrix_size = uintd_to_uint4(matrix_size);
				uint4 _matrix_size_os = uintd_to_uint4(matrix_size_os);
				uint4 _domain_size_grid = uintd_to_uint4(domain_size_grid);

				// Grid oversampled image
				cuFloatComplex *filter_os;
				cudaMalloc((void**)&filter_os, prod(_matrix_size_os)*sizeof(cuFloatComplex));

				compute_deapodization_filter<uint4,float4>( _matrix_size_os, W, 1.0f/W, beta, filter_os );

				// FFT
				if( success )
					success = K2I_ALL( filter_os, _matrix_size_os, 1, false, true );

				// Crop
				if( success )
					crop_image( _matrix_size, _matrix_size_os, out_filterDevPtr, filter_os );
				
				cudaFree(filter_os);

				break;
			}
		default:
			printf("\nCannot calculate deapodization filter of unknown dimension.\n");
			success = false;
			break;
	}

	return success;
}


// Online computation of deapodization filter
template< class UINTd, class FLOATd > __host__ bool
compute_deapodization_filter( UINTd matrix_size_os, float W, float one_over_W, float beta, cuFloatComplex *out_filterDevPtr )
{
	//
	//	Grid fictitious trajectory with a single sample at the origin
	//

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	const unsigned int _blen = 192;
	dim3 dimBlock( _blen, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)prod(matrix_size_os)/(double)_blen), 1, 1 );

	// Invoke kernel
	compute_deapodization_filter_kernel<UINTd, FLOATd><<< dimGrid, dimBlock >>> ( matrix_size_os, W, one_over_W, beta, out_filterDevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'compute_deapodization_filter_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	return true;
}




/*
         Some utility functions
*/


template< class UINTd > void
image_wrap( UINTd matrix_size_os, UINTd matrix_size_wrap, unsigned int number_of_images, bool accumulate, cuFloatComplex *source_image, cuFloatComplex *target_image )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	unsigned int dimBlockX = ((deviceProp.maxThreadsPerBlock/number_of_images)/deviceProp.warpSize)*deviceProp.warpSize;

	// Make sure that gridDim.x becomes integer valued
	while(((prod(matrix_size_os)/get_last_dim(matrix_size_os))%dimBlockX) != 0 ){
		if( dimBlockX<deviceProp.warpSize ){
			printf("\nImplementation error: image_wrap: Cannot reduce the block size below the warp size.\n");
			exit(1);
		}
		else{
			dimBlockX -= deviceProp.warpSize;	
		}
	}

	if(dimBlockX == 0){
		printf("\nImplementation error: image_wrap: Too many coils reduces the block size below the warp size.\n");
		exit(1);
	}

	if( (sizeof(UINTd)!=sizeof(uint2)) && get_last_dim(matrix_size_wrap) != 0 ){
		printf("\nImplementation error: image_wrap: 3D wrapping not yet implemented.\n");
		exit(1);
	}

	dim3 dimBlock( dimBlockX, number_of_images );
	dim3 dimGrid( (prod(matrix_size_os)/get_last_dim(matrix_size_os))/dimBlock.x, get_last_dim(matrix_size_os) ); // No 'ceil'!!!

	// Invoke kernel
	image_wrap_kernel<UINTd><<< dimGrid, dimBlock >>>( matrix_size_os, matrix_size_wrap, number_of_images, accumulate, source_image, target_image );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_wrap_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

__host__ void
cuda_float3_to_float4( float4 *targetDevPtr, float3 *sourceDevPtr, unsigned int number_of_elements )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );
	
	dim3 dimBlock( min((int)number_of_elements, (int)deviceProp.maxThreadsPerBlock), 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)number_of_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	cuda_float3_to_float4_kernel<<< dimGrid, dimBlock >>>( targetDevPtr, sourceDevPtr, number_of_elements );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'cuda_float3_to_float4_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

__host__ void
cuda_uint3_to_uint4( uint4 *targetDevPtr, uint3 *sourceDevPtr, unsigned int number_of_elements )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );
	
	dim3 dimBlock( min((int)number_of_elements, (int)deviceProp.maxThreadsPerBlock), 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)number_of_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	cuda_uint3_to_uint4_kernel<<< dimGrid, dimBlock >>>( targetDevPtr, sourceDevPtr, number_of_elements );
	
	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'cuda_uint3_to_uint4_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


/*
	Kernels
*/

// Kaiser-Bessel convolution kernels

inline __device__ float 
bessi0(float x)
{
	// From numerical recipes in C
	float ax,ans,y;
	if ((ax=fabsf(x)) < 3.75f) 
	{
		y=x/3.75f;
		y*=y;
		ans=1.0f+y*(3.5156229f+y*(3.0899424f+y*(1.2067492f+y*(0.2659732f+y*(0.0360768f+y*0.0045813f)))));
	} 
	else 
	{
		y=3.75f/ax;
		ans=(-0.02057706f+y*(0.02635537f+y*(-0.01647633f+(y*0.00392377f))));
		ans=(__expf(ax)/sqrtf(ax))*(0.39894228f+y*(0.01328592f+y*(0.00225319f+y*(-0.00157565f+y*(0.00916281f+y*ans)))));
	}
	return ans;
}


// Kaiser Bessel according to Beatty et. al. IEEE TMI 2005;24(6):799-808.
// Notice slight difference vs Jackson's formulation, IEEE TMI 1991;10(3):473-478.

inline __device__ float
KaiserBessel( float u, float matrix_size_os_f1, float one_over_W, float beta )
{
  float _tmp = 2.0f*u*one_over_W;
  float tmp = _tmp*_tmp;
  float arg = beta*sqrtf(1.0f-tmp);
  float bessi = bessi0(arg);
  float ret = matrix_size_os_f1*bessi*one_over_W;
  return ret;
}

inline __device__ float
KaiserBessel( float2 u, float2 matrix_size_os, float one_over_W, float beta, bool lastDimFixed )
{
  float phi_x = KaiserBessel( u.x, matrix_size_os.x, one_over_W, beta );
  float phi_y = KaiserBessel( u.y, matrix_size_os.y, one_over_W, beta );

  return phi_x*phi_y;
}

inline __device__ float
KaiserBessel( float3 u, float3 matrix_size_os, float one_over_W, float beta, bool lastDimFixed )
{
  float phi_x = KaiserBessel( u.x, matrix_size_os.x, one_over_W, beta );
  float phi_y = KaiserBessel( u.y, matrix_size_os.y, one_over_W, beta );
  float phi_z;
  if( !lastDimFixed )
	  phi_z = KaiserBessel( u.z, matrix_size_os.z, one_over_W, beta );
  else
	  phi_z = 1.0f;

  return phi_x*phi_y*phi_z;
}

inline __device__ float
KaiserBessel( float4 u, float4 matrix_size_os, float one_over_W, float beta, bool lastDimFixed )
{
  float phi_x = KaiserBessel( u.x, matrix_size_os.x, one_over_W, beta );
  float phi_y = KaiserBessel( u.y, matrix_size_os.y, one_over_W, beta );
  float phi_z = KaiserBessel( u.z, matrix_size_os.z, one_over_W, beta );
  float phi_w;
  if( !lastDimFixed )
	  phi_w = KaiserBessel( u.w, matrix_size_os.w, one_over_W, beta );
  else
	  phi_w = 1.0f;

  return phi_x*phi_y*phi_z*phi_w;
}

// Density compensation

__global__ void 
density_compensate_kernel( unsigned int number_of_samples, cuFloatComplex *out, const cuFloatComplex *in, const float *weights, unsigned int number_of_images )
{
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
	if( index < number_of_samples ){
		const float weight = weights[index];
		for( unsigned int image=0; image<number_of_images; image++ ){		
			const cuFloatComplex tmp = in[image*number_of_samples+index];
			const cuFloatComplex res = make_cuFloatComplex( cuCrealf(tmp)*weight, cuCimagf(tmp)*weight );
			out[image*number_of_samples+index] = res;
		}
	}
}
	

//	Deapodization

// Sizes of image and filter match
__global__ void 
_deapodize_kernel( unsigned int num_elements, cuFloatComplex *imageDevPtr, cuFloatComplex *filterDevPtr, unsigned int number_of_images )
{
	// TODO: use multiplication instead
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){
		for( unsigned int image=0; image<number_of_images; image++ )
			imageDevPtr[image*num_elements+idx] = cuCdivf( imageDevPtr[image*num_elements+idx], filterDevPtr[idx] );
	}
}

// Differently sized image and filter
template< class UINTd > __global__ void 
deapodize_kernel( UINTd image_size, UINTd filter_size, UINTd corner1, UINTd corner2, cuFloatComplex *imageDevPtr, cuFloatComplex *filterDevPtr, unsigned int number_of_images )
{
	// TODO: use multiplication instead
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int num_elements = prod(image_size);
	if( idx < num_elements ){
		const UINTd co = idx_to_co( idx, image_size );
		if( co>=corner1 && co<corner2 ){
			unsigned int filter_idx = co_to_idx(co-corner1, filter_size);
			for( unsigned int image=0; image<number_of_images; image++ )
				imageDevPtr[image*num_elements+idx] = cuCdivf( imageDevPtr[image*num_elements+idx], filterDevPtr[filter_idx] );
		}
	}
}


template< class UINTd, class FLOATd > __global__ void
compute_deapodization_filter_kernel( UINTd matrix_size_os, float W, float one_over_W, float beta, cuFloatComplex *out_DevPtr )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int num_elements = prod(matrix_size_os);
	if( idx < num_elements ){

		// Compute weight from Kaiser-Bessel filter
		const UINTd cell_pos = idx_to_co(idx, matrix_size_os);

		// Sample position ("origin")
		const UINTd sample_pos = matrix_size_os>>1;

		// Calculate the distance between the cell and the sample
		const FLOATd delta = fabsf(uintd_to_floatd(sample_pos)-uintd_to_floatd(cell_pos)); // avoid unsigned subtraction

		// Compute convolution weight. 
		float weight;

		if( weak_greater(delta, W/2.0f ) )

			weight = 0.0f;

		else{ 

			weight = KaiserBessel( delta, uintd_to_floatd(matrix_size_os), one_over_W, beta, false );

			if( !isfinite(weight) )
				weight = 0.0f;
		}

		// Output weight
		out_DevPtr[idx] =  make_cuFloatComplex( weight, 0 );
	}
}


// Image wrap kernels

__inline__ __device__ void 
__wrap_image( uint2 co, uint2 matrix_size_os, uint2 half_wrap, uint2 source_image_size, unsigned int source_image_offset, cuFloatComplex *source_image, cuFloatComplex *out )
{
	// Wrap corners
	if( co.x < half_wrap.x )
		*out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y),source_image_size)+source_image_offset];
	if( co.x >= matrix_size_os.x-half_wrap.x )
		*out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y),source_image_size)+source_image_offset];
	if( co.y < half_wrap.y )
		*out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x, co.y+half_wrap.y+matrix_size_os.y),source_image_size)+source_image_offset];
	if( co.y >= matrix_size_os.y-half_wrap.y )
		*out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x, co.y+half_wrap.y-matrix_size_os.y),source_image_size)+source_image_offset];

	// Wrap edges
	if( co.x < half_wrap.x && co.y < half_wrap.y )
		*out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y+matrix_size_os.y),source_image_size)+source_image_offset];
	if( co.x < half_wrap.x && co.y >= matrix_size_os.y-half_wrap.y )
		*out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y-matrix_size_os.y),source_image_size)+source_image_offset];
	if( co.x >= matrix_size_os.x-half_wrap.x && co.y < half_wrap.y )
		*out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y+matrix_size_os.y),source_image_size)+source_image_offset];
	if( co.x >= matrix_size_os.x-half_wrap.x && co.y >= matrix_size_os.y-half_wrap.y )
		*out += source_image[co_to_idx(make_uint2(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y-matrix_size_os.y),source_image_size)+source_image_offset];
}

__inline__ __device__ void 
__wrap_image( uint3 co, uint3 matrix_size_os, uint3 half_wrap, uint3 source_image_size, unsigned int source_image_offset, cuFloatComplex *source_image, cuFloatComplex *out )
{
	// !!! THIS CODE IS ONLY VALID WHEN THE THERE IS NO CONVOLUTION IN THE 'Z' DIMENSION (k-t SENSE) !!!

	// Wrap corners
	if( co.x < half_wrap.x )
		*out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y, co.z),source_image_size)+source_image_offset];
	if( co.x >= matrix_size_os.x-half_wrap.x )
		*out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y, co.z),source_image_size)+source_image_offset];
	if( co.y < half_wrap.y )
		*out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x, co.y+half_wrap.y+matrix_size_os.y, co.z),source_image_size)+source_image_offset];
	if( co.y >= matrix_size_os.y-half_wrap.y )
		*out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x, co.y+half_wrap.y-matrix_size_os.y, co.z),source_image_size)+source_image_offset];

	// Wrap edges
	if( co.x < half_wrap.x && co.y < half_wrap.y )
		*out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y+matrix_size_os.y, co.z),source_image_size)+source_image_offset];
	if( co.x < half_wrap.x && co.y >= matrix_size_os.y-half_wrap.y )
		*out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x+matrix_size_os.x, co.y+half_wrap.y-matrix_size_os.y, co.z),source_image_size)+source_image_offset];
	if( co.x >= matrix_size_os.x-half_wrap.x && co.y < half_wrap.y )
		*out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y+matrix_size_os.y, co.z),source_image_size)+source_image_offset];
	if( co.x >= matrix_size_os.x-half_wrap.x && co.y >= matrix_size_os.y-half_wrap.y )
		*out += source_image[co_to_idx(make_uint3(co.x+half_wrap.x-matrix_size_os.x, co.y+half_wrap.y-matrix_size_os.y, co.z),source_image_size)+source_image_offset];

}

__inline__ __device__ void 
__wrap_image( uint4 co, uint4 matrix_size_os, uint4 half_wrap, uint4 source_image_size, unsigned int source_image_offset, cuFloatComplex *source_image, cuFloatComplex *out )
{
}

template< class UINTd > __global__ void
image_wrap_kernel( UINTd matrix_size_os, UINTd matrix_size_wrap, unsigned int number_of_images, bool accumulate, cuFloatComplex *source_image, cuFloatComplex *target_image )
{
	const unsigned int image_idx = blockIdx.y*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x;
	const unsigned int image_number = threadIdx.y;
	const unsigned int idx = image_number*gridDim.y*gridDim.x*blockDim.x+image_idx;

	if( image_idx < prod(matrix_size_os) ){

		const UINTd source_image_size = matrix_size_os+matrix_size_wrap;
		const unsigned int source_image_offset = image_number*prod(source_image_size);

		const UINTd half_wrap = matrix_size_wrap>>1;
		const UINTd co = idx_to_co( image_idx, matrix_size_os );

		cuFloatComplex out = source_image[co_to_idx(co+half_wrap,source_image_size)+source_image_offset];

		if( accumulate )
			out += target_image[idx];

		// Wrap image
		__wrap_image( co, matrix_size_os, half_wrap, source_image_size, source_image_offset, source_image, &out ); 
		
		target_image[idx] = out;
	}
}

// Private building block kernels


__global__ void
cuda_float3_to_float4_kernel( float4 *targetDevPtr, float3 *sourceDevPtr, unsigned int number_of_elements )
{
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
	if( index < number_of_elements ){

		const float3 in = sourceDevPtr[index];
		float4 out;

		out.x = in.x;
		out.y = in.y;
		out.z = in.z;
		out.w = 0.0f;

		targetDevPtr[index] = out;
	}
}

__global__ void
cuda_uint3_to_uint4_kernel( uint4 *targetDevPtr, uint3 *sourceDevPtr, unsigned int number_of_elements )
{
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
	if( index < number_of_elements ){

		const uint3 in = sourceDevPtr[index];
		uint4 out;

		out.x = in.x;
		out.y = in.y;
		out.z = in.z;
		out.w = 0.0f;

		targetDevPtr[index] = out;
	}
}


/*
     Kernels common to the different NFFT and NFFT_H scenarios
     The actual NFFT and NFFT_H kernels are in individual files (NFFT_NFFT*.cu).
*/


// Return UINTd strips direction for thread 'globalThreadId'

template <class UINTd, bool CORE> inline __device__ UINTd
NFFT_get_strip_dir( unsigned int globalThreadId )
{

	// A good compiler/optimizer (presumably 'nvcc') will detect that this switch can be fully "unrolled" at compile time

	UINTd dir;

	switch( sizeof(UINTd) )
	{
	case sizeof(uint2):
		{
			uint2 tmp = tex1Dfetch( TEX_STRIPS_DIR_UI2_NFFT(CORE), globalThreadId );
			dir = *((UINTd*) &tmp);
		}
		break;
	case sizeof(uint3):
		{
			// uint3 textures doesn't exist. We use uint4 for the fetch instead.
			uint4 tmp = tex1Dfetch( TEX_STRIPS_DIR_UI4_NFFT(CORE), globalThreadId );
			dir = *((UINTd*) &tmp);
		}
		break;
	case sizeof(uint4):
		{
			uint4 tmp = tex1Dfetch( TEX_STRIPS_DIR_UI4_NFFT(CORE), globalThreadId );
			dir = *((UINTd*) &tmp);
		}
		break;
	}

	return dir;
}


// Return UINTd strip origin for strip 'stripIdx'

template <class UINTd, bool CORE> inline __device__ UINTd
NFFT_get_strip_origin( unsigned int stripIdx )
{

	// A good compiler/optimizer (presumably 'nvcc') will detect that this switch can be fully "unrolled" at compile time

	UINTd origin;

	switch( sizeof(UINTd) )
	{
	case sizeof(uint2):
		{
			uint2 tmp = tex1Dfetch( TEX_STRIPS_ORIGINS_UI2_NFFT(CORE), stripIdx );						
			origin = *((UINTd*) &tmp);
		}
		break;
	case sizeof(uint3):
		{
			// uint3 textures doesn't exist. We use uint4 for the fetch instead.
			uint4 tmp = tex1Dfetch( TEX_STRIPS_ORIGINS_UI4_NFFT(CORE), stripIdx );						
			origin = *((UINTd*) &tmp);
		}
		break;
	case sizeof(uint4):
		{
			uint4 tmp = tex1Dfetch( TEX_STRIPS_ORIGINS_UI4_NFFT(CORE), stripIdx );						
			origin = *((UINTd*) &tmp);
		}
		break;
	}

	return origin;
}


// Return FLOATd sample position of sample index 'sampleIdx'

template <class FLOATd, bool CORE> inline __device__ FLOATd
NFFT_get_sample_position( unsigned int sampleIdx )
{

	// A good compiler/optimizer (presumably 'nvcc') will detect that this switch can be fully "unrolled" at compile time

	FLOATd sample_position;

	switch( sizeof(FLOATd) )
	{
	case sizeof(float2):
		{
			float2 tmp = tex1Dfetch( TEX_SAMPLE_POSITIONS_F2_NFFT(CORE), sampleIdx );
			sample_position = *((FLOATd*) &tmp);
		}
		break;
	case sizeof(float3):
		{
			// float3 textures doesn't exist. We use float4 for the fetch instead.
			float4 tmp = tex1Dfetch( TEX_SAMPLE_POSITIONS_F4_NFFT(CORE), sampleIdx );
			sample_position = *((FLOATd*) &tmp);
		}
		break;
	case sizeof(float4):
		{
			float4 tmp = tex1Dfetch( TEX_SAMPLE_POSITIONS_F4_NFFT(CORE), sampleIdx );
			sample_position = *((FLOATd*) &tmp);
		}
		break;
	}
	return sample_position;
}

template <class FLOATd, bool CORE> inline __device__ FLOATd
NFFT_H_get_sample_position( unsigned int sampleIdx )
{

	// A good compiler/optimizer (presumably 'nvcc') will detect that this switch can be fully "unrolled" at compile time

	FLOATd pos;

	switch( sizeof(FLOATd) )
	{
	case sizeof(float2):
		{
			float2 tmp = tex1Dfetch( TEX_SAMPLE_POSITIONS_F2_NFFT_H(CORE), sampleIdx );
			pos = *((FLOATd*) &tmp);
		}
		break;
	case sizeof(float3):
		{
			// float3 textures doesn't exist. We use float4 for the fetch instead.
			float4 tmp = tex1Dfetch( TEX_SAMPLE_POSITIONS_F4_NFFT_H(CORE), sampleIdx );
			pos = *((FLOATd*) &tmp);
		}
		break;
	case sizeof(float4):
		{
			float4 tmp = tex1Dfetch( TEX_SAMPLE_POSITIONS_F4_NFFT_H(CORE), sampleIdx );
			pos = *((FLOATd*) &tmp);
		}
		break;
	}

	return pos;
}

 inline __device__ void float2_to_floatd_with_temporal_filling( float2 *res, float2 ui, float time )
{
  (*res).x = ui.x;
  (*res).y = ui.y;
}

 inline __device__ void float2_to_floatd_with_temporal_filling( float3 *res, float2 ui, float time )
{
  (*res).x = ui.x;
  (*res).y = ui.y;
  (*res).z = time;
}

inline __device__ void float2_to_floatd_with_temporal_filling( float4 *res, float2 ui, float time )
{
  (*res).x = ui.x;
  (*res).y = ui.y;
  (*res).z = time;
  (*res).w = 0.0f;
}

// Local version for the NFFT and included .cu files
template <class UINTd, class FLOATd, char TYPE> __inline__ __device__ FLOATd
cuda_compute_radial_sample_position( unsigned int _sampleIdx, float alpha, UINTd bias, UINTd bias_os )
{
	FLOATd sample_pos;
	float cos_angle, sin_angle;

	switch(TYPE){

		// We use the "fast math" version of sine and cosine. 
		// Otherwise the compiler resolves to local memory (nvcc Cuda 2.0)

	case 1:

		// Golden angle
		{
			const float angle_step = CUDART_PI_F/((sqrtf(5.0f) + 1.0f)/2.0f);
			const float sampleIdx = uint2float(_sampleIdx);
			const float frame = floorf((sampleIdx+0.5f)*__one_over_num_samples_per_frame);
			const float e1 = floorf((sampleIdx+0.5f)*__one_over_num_samples_per_projection);
			__sincosf( (e1+__angular_offset)*angle_step*__gc_factor, &sin_angle, &cos_angle );
			const float2 _sample_pos = make_float2( (((sampleIdx-e1*__num_samples_per_projection)*__one_over_radial_oversampling)-bias.x)*cos_angle*alpha+bias_os.x, (((sampleIdx-e1*__num_samples_per_projection)*__one_over_radial_oversampling)-bias.y)*sin_angle*alpha+bias_os.y );
			float2_to_floatd_with_temporal_filling( &sample_pos, _sample_pos, frame);
		}

		break;

	case 2:

		// Fixed angle radial projections
		{
			const float sampleIdx = uint2float(_sampleIdx);
			const float frame = floorf((sampleIdx+0.5f)*__one_over_num_samples_per_frame);
			const float local_sample_idx = sampleIdx - frame*__num_samples_per_frame;
			const float e1 = floorf((local_sample_idx+0.5f)*__one_over_num_samples_per_projection);
			const float rotation = fmodf(fmodf(frame+__angular_offset, __frames_per_rotation_cycle)*__rotation_gap, __frames_per_rotation_cycle);
			__sincosf( e1*CUDART_PI_F*__one_over_num_projections_per_frame+rotation*__interframe_rotation, &sin_angle, &cos_angle );
			const float2 _sample_pos = make_float2( (((local_sample_idx-e1*__num_samples_per_projection)*__one_over_radial_oversampling)-bias.x)*cos_angle*alpha+bias_os.x, (((local_sample_idx-e1*__num_samples_per_projection)*__one_over_radial_oversampling)-bias.y)*sin_angle*alpha+bias_os.y );
			float2_to_floatd_with_temporal_filling( &sample_pos, _sample_pos, frame);
		}

		break;

	default:
		break;
	}

	return sample_pos;	
}

template <class UINTd, class FLOATd, char TYPE> __global__ void 
compute_radial_sample_positions_kernel( UINTd bias, UINTd bias_os, FLOATd matrix_size_os_f, FLOATd bias_os_f, unsigned int number_of_samples, float alpha, FLOATd *co, bool remove_os, bool normalize )
{
	const unsigned int index = blockIdx.x*blockDim.x + threadIdx.x;
	
	if( index < number_of_samples ){

		 FLOATd pos = cuda_compute_radial_sample_position<UINTd,FLOATd,TYPE>( index, alpha, bias, bias_os );
		 if( remove_os ){
			 if( normalize )
				co[index] = (pos-bias_os_f)/matrix_size_os_f;
			else
				co[index] = pos/alpha;
		 }
		else
			// Ignoring normalization option in this case
			co[index] = pos;
	}
}

template <class UINTd, class FLOATd, char TYPE> __host__
FLOATd* compute_radial_sample_positions( unsigned int number_of_samples, UINTd matrix_size, UINTd matrix_size_os, float alpha, unsigned int samples_per_projection, unsigned int projections_per_frame, unsigned int angular_offset, unsigned int frames_per_rotation, float gc_factor, bool remove_os, bool normalize )
{
	// Find dimensions of grid/blocks.
	const unsigned int block_size = min(256,number_of_samples);
	dim3 dimBlock( block_size, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)number_of_samples/dimBlock.x), 1, 1 );

	// Allocate space for result
	FLOATd *co;
	cudaMalloc((void**)&co, number_of_samples*sizeof(FLOATd));
	
	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'compute_radial_sample_positions' 1: %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}

	// Invoke kernel
	_nfft_NFFT_set_constants<TYPE>( matrix_size.x, samples_per_projection, projections_per_frame, angular_offset, frames_per_rotation, gc_factor );
	compute_radial_sample_positions_kernel<UINTd,FLOATd,TYPE><<< dimGrid, dimBlock >>> ( matrix_size>>1, matrix_size_os>>1, cast_uintd_to_floatd(matrix_size_os),cast_uintd_to_floatd(matrix_size_os>>1), number_of_samples, alpha, co, remove_os, normalize );

	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'compute_radial_sample_positions_kernel' 2: %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}

	return co;
}

template <class UINTd, class FLOATd, char TYPE> __host__
FLOATd* compute_radial_sample_positions( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE> *plan, bool remove_os, bool normalize )
{
	return compute_radial_sample_positions<UINTd,FLOATd,TYPE>( plan->number_of_samples, plan->matrix_size, plan->matrix_size_os, plan->alpha, plan->samples_per_projection, plan->projections_per_frame, plan->angular_offset, plan->frames_per_rotation, plan->gc_factor, remove_os, normalize );
}


template <class UINTd, class FLOATd> __host__
float* estimate_dcw( mr_recon::NFFT_iteration_plan<UINTd, FLOATd, 0> *plan, FLOATd *trajDevPtr )
{
	// Allocate device memory
	cuFloatComplex *samples, *image;
	float *dcw;
	cudaMalloc((void**)&samples, plan->number_of_samples*sizeof(cuFloatComplex));
	cudaMalloc((void**)&image, prod(plan->matrix_size_os)*sizeof(cuFloatComplex));
	cudaMalloc((void**)&dcw, plan->number_of_samples*sizeof(float));

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'estimate_dcw': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}

	// Initialize samples to complex(1,0)
	clear_image( plan->number_of_samples, make_cuFloatComplex(1.0f, 0.0f), samples );

	// Convolve onto grid
	bool success = 
		NFFT_iteration_convolve_to_image( plan, samples, image );

	// "normalize" - limit max value 
	normalize_max_length( prod(plan->matrix_size_os), image );

	// Convolve onto original sample positions
	if( success ) success =
		NFFT_iteration_convolve_to_samples( plan, samples, image );

	// "normalize" - limit max value 
	normalize_max_length( plan->number_of_samples, samples );

	if( !success ){
		cudaFree(samples);
		cudaFree(image);
		cudaFree(dcw);
		return 0x0;
	}

	// Compute image magnitudes. 
	image_modulus( samples, dcw, plan->number_of_samples, false );

	cudaFree(samples);
	cudaFree(image);

	return dcw;
}


/*
    --------------------
    Code modularization.
	--------------------
	
	Let domain_size_spatial denote the number of convolution target points per coil for each processor.
	Let domain_size_coils denote the number of coils the convolution should be applied to.
	This corresponds to the domain_size_* variables passed to preprocess_NFFT (preprocess.hpp). 

	For performance reasons we have dedicated code paths for each of several scenarios:

	  - (domain_size_spatial > 1) and (domain_size_coils > 1):   GENERIC    code path. 
	  - (domain_size_spatial = 1) and (domain_size_coils >= 1):  PARALLEL   code path. 
	  - (domain_size_spatial > 1) and (domain_size_coils = 1):   SEQUENTIAL code path. 
       
	Furthermore, the generic path is only semigeneric if the spatial variation occurs only in one dimension.
	
	To make our code more readable, we have forked out our implementation of each path in individual files.
	Each combination is #include'd below.

	The path is determined in the NFFT_convolve_to_image and NFFT_convolve_to_samples functions below.
*/	

const unsigned int my_warpSize = 32;

//#include "NFFT_NFFT_generic.cu"
#include "NFFT_NFFT_parallel.cu"
//#include "NFFT_NFFT_sequential.cu"

//#include "NFFT_NFFT_H_generic.cu"
//#include "NFFT_NFFT_H_semigeneric.cu"
#include "NFFT_NFFT_H_parallel.cu"
//#include "NFFT_NFFT_H_sequential.cu"


/*
	Template instantion
*/

#include "NFFT_template_instantiation.cu"
