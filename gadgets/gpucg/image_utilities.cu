#include "NFFT.hcu"
#include "image_utilities.hcu"
#include "int_util.hcu"
#include "uint_util.hcu"

#include <stdio.h>
#include <cublas.h>
#include <cuComplex.h>

texture<float2, 1> permute_tex;

/*
	Image modulus
*/

__global__ void 
image_modulus_kernel( cuFloatComplex *devPtr_in, float *devPtr_out, unsigned int number_of_elements )
{
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx<number_of_elements ){
	  cuFloatComplex val = devPtr_in[idx];    
	  devPtr_out[idx] = sqrtf(val.x*val.x+val.y*val.y);
  }
}

__host__ void
image_modulus( cuFloatComplex *devPtr_in, float *devPtr_out, unsigned int number_of_elements, bool normalize, float normalize_scale )
{
	dim3 blockDim(512,1,1);
	dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x), 1, 1 );

	// Make modulus image
	image_modulus_kernel<<< gridDim, blockDim >>>( devPtr_in, devPtr_out, number_of_elements );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_modulus_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	if( normalize ){

		/*
			Scale using cublas. 
			Remember cublasIsamax uses 1-based indexing!!!
		*/

		//Find the maximum value in the array
		int max_idx = cublasIsamax (number_of_elements, devPtr_out, 1);

		//Copy that value back to host memory
		float max_val;
		cudaMemcpy(&max_val, (devPtr_out+max_idx-1), sizeof(float), cudaMemcpyDeviceToHost);

		//Scale the array (2.5 is an "arbitrary" scaling constant)
		cublasSscal( number_of_elements, normalize_scale/max_val, devPtr_out, 1 );
	}
}


/*
	Image normalize
*/

__host__ void 
image_normalize( unsigned int number_of_elements, float *devPtr )
{
	/*
		Scale using cublas. 
		Remember cublasIsamax uses 1-based indexing!!!
	*/

	//Find the maximum value in the array
	int max_idx = cublasIsamax (number_of_elements, devPtr, 1);

	//Copy that value back to host memory
	float max_val;
	cudaMemcpy(&max_val, (devPtr+max_idx-1), sizeof(float), cudaMemcpyDeviceToHost);

	// printf("\nMax index/val: %d/%f\n", max_idx, max_val );

	//Scale the array (2.5 is an "arbitrary" scaling constant)
	cublasSscal( number_of_elements, 2.5f/max_val, devPtr, 1 );
}

__host__ void 
normalize_max_length( unsigned int number_of_elements, cuFloatComplex *devPtr )
{
	// Base scaling on magnitudes image (largest element)

	float *tmp;
	cudaMalloc((void**)&tmp, number_of_elements*sizeof(cuFloatComplex));

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_set_max_length': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout); 
		exit(1);
	}

	image_modulus( devPtr, tmp, number_of_elements, false );

	/*
		Scale using cublas. 
		Remember cublasIsamax uses 1-based indexing!!!
	*/

	//Find the maximum value in the array
	int max_idx = cublasIsamax (number_of_elements, tmp, 1);

	//Copy that value back to host memory
	float max_val;
	cudaMemcpy(&max_val, (tmp+max_idx-1), sizeof(float), cudaMemcpyDeviceToHost);

	// printf("\nMax index/val: %d/%f\n", max_idx, max_val );

	//Scale the array (2.5 is an "arbitrary" scaling constant)
//DOESN'T WORK!	cublasCscal( number_of_elements, make_cuFloatComplex(1.0f/max_val,1.0f/max_val), devPtr, 1 );
	complex_image_scale( 1.0f/max_val, number_of_elements, devPtr );
	
	cudaFree(tmp);
}


/* 
	Clear image
*/

template< class T > __global__ void 
clear_image_kernel( unsigned int num_elements, T val, T *imageDevPtr )
{
	const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){
		imageDevPtr[idx] = val;
	}
}

template< class T > __host__ void 
clear_image( unsigned int num_elements, T val, T *imageDevPtr )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	int gridX, gridY;
	int total_blocks = (int) ceil((double)num_elements/deviceProp.maxThreadsPerBlock);

	if(total_blocks<=deviceProp.maxGridSize[0])
	{
		gridX = total_blocks;
		gridY = 1;
	}
	else{
		gridX = total_blocks;
		gridY = 1;
		while( gridX>gridY ){
			gridX=gridX>>1;
			gridY=gridY<<1;
		}
		if(gridX*gridY<total_blocks){
			gridX=(total_blocks/gridY)+((total_blocks%gridY)?1:0);
		}
	}

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( gridX, gridY, 1 );

	// Invoke kernel
	clear_image_kernel<T><<< dimGrid, dimBlock >>>( num_elements, val, imageDevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'clear_image_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


/*
	Border fill image (square)
*/

template < class T, class UINTd> __global__ void 
cuda_border_fill_image_kernel( UINTd matrix_size, UINTd matrix_size_os, T value, T *imageDevPtr, unsigned int number_of_images )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int num_elements = prod(matrix_size_os);

	if( idx < num_elements ){
		const UINTd co = idx_to_co( idx, matrix_size_os );
		const UINTd corner1 = (matrix_size_os-matrix_size)>>1;
		const UINTd corner2 = matrix_size_os-corner1;
		if( weak_less(co,corner1) || weak_greater_equal(co,corner2) )
			for( unsigned int image=0; image<number_of_images; image++ )
				imageDevPtr[image*num_elements+idx] = value;
	}
}

template< class T, class UINTd > __host__ void 
cuda_border_fill_image( UINTd matrix_size, UINTd matrix_size_os, T value, T *imageDevPtr, unsigned int number_of_images )
{
	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)prod(matrix_size_os)/dimBlock.x), 1, 1 );

	// Invoke kernel
	cuda_border_fill_image_kernel<T,UINTd><<< dimGrid, dimBlock >>> ( matrix_size, matrix_size_os, value, imageDevPtr, number_of_images );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'cuda_border_fill_image_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


/*
	Border fill image (circular)
*/
__device__ void 
__arrange( uint2 matrix_size, uint2 &co )
{
	if( co.x < (matrix_size.x>>1) ) 
		co.x = matrix_size.x-co.x;
	if( co.y < (matrix_size.y>>1) ) 
		co.y = matrix_size.y-co.y;
}

__device__ void 
__arrange( uint3 matrix_size, uint3 &co )
{
	if( co.x < (matrix_size.x>>1) ) 
		co.x = matrix_size.x-co.x;
	if( co.y < (matrix_size.y>>1) ) 
		co.y = matrix_size.y-co.y;
	if( co.z < (matrix_size.z>>1) ) 
		co.z = matrix_size.z-co.z;
}

__device__ void 
__arrange( uint4 matrix_size, uint4 &co )
{
	if( co.x < (matrix_size.x>>1) ) 
		co.x = matrix_size.x-co.x;
	if( co.y < (matrix_size.y>>1) ) 
		co.y = matrix_size.y-co.y;
	if( co.z < (matrix_size.z>>1) ) 
		co.z = matrix_size.z-co.z;
	if( co.w < (matrix_size.w>>1) ) 
		co.w = matrix_size.w-co.w;
}

template < class T, class UINTd > __global__ void 
cuda_border_fill_image_kernel( UINTd matrix_size, float radius, T value, T *imageDevPtr, unsigned int number_of_images )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int num_elements = prod(matrix_size);

	if( idx < num_elements ){
		UINTd co = idx_to_co( idx, matrix_size );
		__arrange( matrix_size, co ); // Make sure that co-(matrix_size>>1) is positive
		const UINTd tmp = co-(matrix_size>>1);
		const unsigned int sq_dist = dot(tmp,tmp);
		if( sq_dist>radius*radius ) 
			for( unsigned int image=0; image<number_of_images; image++ )
				imageDevPtr[image*num_elements+idx] = value;
	}
}

template< class T, class UINTd > __host__ void 
cuda_border_fill_image( UINTd matrix_size_os, float radius, T value, T *imageDevPtr, unsigned int number_of_images )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)prod(matrix_size_os)/dimBlock.x), 1, 1 );

	// Invoke kernel
	cuda_border_fill_image_kernel<T,UINTd><<< dimGrid, dimBlock >>> ( matrix_size_os, radius, value, imageDevPtr, number_of_images );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'cuda_border_fill_image_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


/*
	Image crop
*/

template< class T, class UINTd > __global__ void
crop_image_kernel( UINTd out_matrix_size, UINTd in_matrix_size, T *out_imageDevPtr, T *in_imageDevPtr, unsigned int number_of_images )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int num_elements = prod(out_matrix_size);

	if( idx < num_elements ){
		const UINTd co = idx_to_co( idx, out_matrix_size );
		const UINTd co_os = co + ((in_matrix_size-out_matrix_size)>>1);
		const unsigned int source_idx = co_to_idx(co_os, in_matrix_size);
		const unsigned int source_elements = prod(in_matrix_size);
		for( unsigned int image=0; image<number_of_images; image++ )
			out_imageDevPtr[image*num_elements+idx] = in_imageDevPtr[image*source_elements+source_idx];
	}
}

template< class T, class UINTd > 
__host__ void crop_image( UINTd out_matrix_size, UINTd in_matrix_size, T *out_imageDevPtr, T *in_imageDevPtr, unsigned int number_of_images )
{
	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)prod(out_matrix_size)/(double)deviceProp.maxThreadsPerBlock), 1, 1 );

	// Invoke kernel
	crop_image_kernel<T,UINTd><<< dimGrid, dimBlock >>> ( out_matrix_size, in_matrix_size, out_imageDevPtr, in_imageDevPtr, number_of_images );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'crop_image_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

inline __device__ float
gaussian( float squared_u, float minus_half_over_sigmaSquared )
{
	return expf( squared_u * minus_half_over_sigmaSquared );
}

inline __device__ void __image_multiply( float *ptr, float val )
{
	*ptr *= val;
}

inline __device__ void __image_multiply( cuFloatComplex *ptr, float val )
{
	cuFloatComplex tmp = make_cuFloatComplex( cuCrealf(*ptr)*val, cuCimagf(*ptr)*val );
	*ptr = tmp;
}

inline __device__ void __image_multiply( cuFloatComplex *ptr, cuFloatComplex val )
{
	*ptr = cuCmulf( *ptr, val );
}

// Multiply image with Gaussian
template< class T, class UINTd >
__global__ void image_multiply_gaussian_kernel( UINTd matrix_size, UINTd center, float minus_half_over_sigmaSquared, T *imageDevPtr, unsigned int number_of_images = 1 )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int num_elements = prod(matrix_size);

	if( idx < num_elements ){
		const UINTd co = idx_to_co( idx, matrix_size );
		const UINTd diff = co-center;
		float _gaussian = gaussian( dot(diff,diff), minus_half_over_sigmaSquared );

		for( unsigned int image=0; image<number_of_images; image++ )
			__image_multiply( &imageDevPtr[image*num_elements+idx], _gaussian );
	}
}

template< class T, class UINTd >
__host__ void image_multiply_gaussian( UINTd matrix_size, float sigma, T *imageDevPtr, unsigned int number_of_images )
{
	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)prod(matrix_size)/(double)dimBlock.x), 1, 1 );

	// Invoke kernel
	image_multiply_gaussian_kernel<T,UINTd><<< dimGrid, dimBlock >>>( matrix_size, matrix_size>>1, -1.0f/(2.0f*sigma*sigma), imageDevPtr, number_of_images );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_multiply_gaussian_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


/*
	Add two images
*/

inline __device__ float __image_add( float val1, float val2 )
{
	return val1 + val2;
}

inline __device__ cuFloatComplex __image_add( cuFloatComplex val1, cuFloatComplex val2 )
{
	return cuCaddf( val1, val2 );
}

inline __device__ float __image_sub( float val1, float val2 )
{
	return val1 - val2;
}

inline __device__ cuFloatComplex __image_sub( cuFloatComplex val1, cuFloatComplex val2 )
{
	return cuCsubf( val1, val2 );
}

template< class T > __global__ void 
add_images_kernel( unsigned int num_elements, T *target, T *source1, T *source2  )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements )
		target[idx] = __image_add( source1[idx], source2[idx] );
}

template< class T > __host__ void 
add_images( unsigned int num_elements, T *targetDevPtr, T *source1DevPtr,  T *source2DevPtr  )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	add_images_kernel<T><<< dimGrid, dimBlock >>>( num_elements, targetDevPtr, source1DevPtr, source2DevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'add_images_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

/*
	Multiply two images
*/

inline __device__ float __image_mult( float val1, float val2 )
{
	return val1 * val2;
}

inline __device__ cuFloatComplex __image_mult( cuFloatComplex val1, cuFloatComplex val2 )
{
	return cuCmulf( val1, val2 );
}


template< class T > __global__ void 
multiply_images_kernel( unsigned int num_elements, T *target, T *source1, T *source2  )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements )
		target[idx] = __image_mult( source1[idx], source2[idx] );
}

template< class T > __host__ void 
multiply_images( unsigned int num_elements, T *targetDevPtr, T *source1DevPtr,  T *source2DevPtr  )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	multiply_images_kernel<T><<< dimGrid, dimBlock >>>( num_elements, targetDevPtr, source1DevPtr, source2DevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'multiply_images_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

inline __device__ void __zero_element( float *ptr )
{
	*ptr = 0.0f;
}

inline __device__ void __zero_element( cuFloatComplex *ptr )
{
	*ptr = make_cuFloatComplex( 0.0f, 0.0f );
}

template< class T > __global__ void 
add_images_kernel( unsigned int num_elements, T *target, T *source, unsigned int num_images  )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){

		T out;

		__zero_element( &out );

		for( unsigned int i=0; i<num_images; i++ )
			out = __image_add( out, source[i*num_elements+idx] );

		target[idx] = out;
	}
}

template< class T > void 
add_images( unsigned int num_elements, T *targetDevPtr, T *sourceDevPtr, unsigned int num_source_images )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	add_images_kernel<T><<< dimGrid, dimBlock >>>( num_elements, targetDevPtr, sourceDevPtr, num_source_images );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'add_images_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

template< class T > __global__ void 
subtract_images_kernel( unsigned int num_elements, T *target, T *source1, T *source2  )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements )
		target[idx] = __image_sub( source1[idx], source2[idx] );
}

template< class T > __host__ void 
subtract_images( unsigned int num_elements, T *targetDevPtr, T *source1DevPtr,  T *source2DevPtr  )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	subtract_images_kernel<T><<< dimGrid, dimBlock >>>( num_elements, targetDevPtr, source1DevPtr, source2DevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'subtract_images_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


/* 
	Square image
*/

template< class T > __global__ void 
square_image_kernel( unsigned int num_elements, T *imageDevPtr )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if( idx < num_elements ){
		T source = imageDevPtr[idx];
		T out = source;

		__image_multiply( &out, source );

		imageDevPtr[idx] = out;
	}
}

template< class T > __host__ void 
square_image( unsigned int num_elements, T *imageDevPtr )
{
	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/(double)dimBlock.x), 1, 1 );

	// Invoke kernel
	square_image_kernel<T><<< dimGrid, dimBlock >>> ( num_elements, imageDevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'square_image_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

/* 
	Squareroot of image
*/

__global__ void 
squareroot_modulus_image_kernel( unsigned int num_elements, float *imageDevPtr )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if( idx < num_elements )
		imageDevPtr[idx] = sqrtf(imageDevPtr[idx]);
}

__host__ void 
squareroot_modulus_image( unsigned int num_elements, float *imageDevPtr )
{
	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/(double)dimBlock.x), 1, 1 );

	// Invoke kernel
	squareroot_modulus_image_kernel<<< dimGrid, dimBlock >>> ( num_elements, imageDevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'squareroot_modulus_image_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

// Reciprocal image
__global__ void 
image_reciprocal_kernel( unsigned int num_elements, float *imageDevPtr )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){
		imageDevPtr[idx] = 1.0f/imageDevPtr[idx];
	}
}

// Reciprocal image
__host__ void 
image_reciprocal( unsigned int num_elements, float *imageDevPtr )
{

	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;  
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	image_reciprocal_kernel<<< dimGrid, dimBlock >>>( num_elements, imageDevPtr );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_reciprocal_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

template< class UINTd > __global__ void
image_copy_with_zero_fill_kernel( UINTd matrix_size, UINTd matrix_size_os, cuFloatComplex *image_DevPtr, cuFloatComplex *image_os_DevPtr, unsigned int number_of_images )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int num_elements = prod(matrix_size_os);
	if( idx < num_elements ){
		const UINTd co_os = idx_to_co( idx, matrix_size_os );
		const UINTd corner1 = (matrix_size_os-matrix_size)>>1;
		const UINTd corner2 = matrix_size_os-corner1;
		if( weak_less(co_os,corner1) || weak_greater_equal(co_os,corner2) ){
			for( unsigned int image=0; image<number_of_images; image++ )
				image_os_DevPtr[image*num_elements+idx] = make_cuFloatComplex( 0.0f, 0.0f );
		}
		else{
			const unsigned int source_elements = prod(matrix_size);
			const unsigned int source_idx = co_to_idx(co_os-corner1, matrix_size);
			for( unsigned int image=0; image<number_of_images; image++ )
				image_os_DevPtr[image*num_elements+idx] = image_DevPtr[image*source_elements+source_idx];
		}
	}
}

template< class UINTd > __host__ void
image_copy_with_zero_fill( UINTd matrix_size, UINTd matrix_size_os, cuFloatComplex *image_DevPtr, cuFloatComplex *image_os_DevPtr, unsigned int number_of_images )
{
	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)prod(matrix_size_os)/dimBlock.x), 1, 1 );

	// Invoke kernel
	image_copy_with_zero_fill_kernel<UINTd><<< dimGrid, dimBlock >>> ( matrix_size, matrix_size_os, image_DevPtr, image_os_DevPtr, number_of_images );	

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_copy_with_zero_fill_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

__global__ void 
complex2real_kernel( unsigned int num_elements, cuFloatComplex *_in, float *out_real, float *out_imag )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < num_elements ){
		cuFloatComplex in = _in[idx];
		out_real[idx] = cuCrealf(in);
		out_imag[idx] = cuCimagf(in);
	}
}

// Split complex image into two real images
__host__ void 
complex2real( unsigned int num_elements, cuFloatComplex *in, float *out_real, float *out_imag )
{
	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)num_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	complex2real_kernel<<< dimGrid, dimBlock >>> ( num_elements, in, out_real, out_imag );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_copy_with_zero_fill_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

// Complex scale
__global__ void 
complex_image_scale_kernel( float scale, unsigned int number_of_elements, cuFloatComplex *image )
{
	const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if( idx < number_of_elements ){
		cuFloatComplex in = image[idx];
		image[idx] = make_cuFloatComplex( scale*cuCrealf(in), scale*cuCimagf(in));
	}
}

__host__ void 
complex_image_scale( float scale, unsigned int number_of_elements, cuFloatComplex *image )
{
	// Find dimensions of grid/blocks.

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties( &deviceProp, _convolution_device );

	dim3 dimBlock( deviceProp.maxThreadsPerBlock, 1, 1 );
	dim3 dimGrid( (unsigned int) ceil((double)number_of_elements/dimBlock.x), 1, 1 );

	// Invoke kernel
	complex_image_scale_kernel<<< dimGrid, dimBlock >>> ( scale, number_of_elements, image );

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'complex_image_scale': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

__global__ void image_permute_kernel( cuFloatComplex* data_in, cuFloatComplex* data_out, uint4 dim, int num_shifts )
{
	//This is the current pixel number
	unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	// Do nothing if index is out of range
	if( idx < prod(dim) ){
		uint4 new_dim = shift( dim, num_shifts );
		uint4 co = idx_to_co(idx, new_dim);
		co = shift( co, -num_shifts );
        data_out[idx] = tex1Dfetch( permute_tex, co_to_idx(co,dim) );
	}
}

void image_permute(cuFloatComplex* data_in, cuFloatComplex* data_out, uint4 dim, int num_shifts)
{
    dim3 dimBlock(512,1,1);
    unsigned int numCellsInGrid = ceil( ((double)prod(dim))/dimBlock.x );
    dim3 dimGrid(numCellsInGrid, 1, 1);

    // setup texture
    const cudaChannelFormatDesc desc_float2 = cudaCreateChannelDesc<float2>();
    cudaBindTexture( 0, &permute_tex, (float2*)data_in, &desc_float2, prod(dim)*sizeof(float2));
	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_permute on txture binding': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

	// Invoke kernel
	image_permute_kernel<<< dimGrid, dimBlock >>> ( data_in, data_out, dim, num_shifts );
	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_permute on kernel launch': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}

    cudaUnbindTexture( permute_tex );
	err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'image_permute on unbind': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}

// Instanciation

template void cuda_border_fill_image( uint2, uint2, cuFloatComplex, cuFloatComplex*, unsigned int );
template void cuda_border_fill_image( uint3, uint3, cuFloatComplex, cuFloatComplex*, unsigned int );
template void cuda_border_fill_image( uint4, uint4, cuFloatComplex, cuFloatComplex*, unsigned int );
template void cuda_border_fill_image( uint2, uint2, float, float*, unsigned int );
template void cuda_border_fill_image( uint3, uint3, float, float*, unsigned int );
template void cuda_border_fill_image( uint4, uint4, float, float*, unsigned int );

template void cuda_border_fill_image<cuFloatComplex, uint2>( uint2, float, cuFloatComplex, cuFloatComplex*, unsigned int );
template void cuda_border_fill_image<cuFloatComplex, uint3>( uint3, float, cuFloatComplex, cuFloatComplex*, unsigned int );
template void cuda_border_fill_image<cuFloatComplex, uint4>( uint4, float, cuFloatComplex, cuFloatComplex*, unsigned int );
template void cuda_border_fill_image<float, uint2>( uint2, float, float, float*, unsigned int );
template void cuda_border_fill_image<float, uint3>( uint3, float, float, float*, unsigned int );
template void cuda_border_fill_image<float, uint4>( uint4, float, float, float*, unsigned int );

template void crop_image( uint2, uint2, cuFloatComplex*, cuFloatComplex*, unsigned int );
template void crop_image( uint3, uint3, cuFloatComplex*, cuFloatComplex*, unsigned int );
template void crop_image( uint4, uint4, cuFloatComplex*, cuFloatComplex*, unsigned int );
template void crop_image( uint2, uint2, float*, float*, unsigned int );
template void crop_image( uint3, uint3, float*, float*, unsigned int );
template void crop_image( uint4, uint4, float*, float*, unsigned int );

template void image_multiply_gaussian( uint2, float, float*, unsigned int );
template void image_multiply_gaussian( uint3, float, float*, unsigned int );
template void image_multiply_gaussian( uint4, float, float*, unsigned int );

template void image_multiply_gaussian( uint2, float, cuFloatComplex*, unsigned int );
template void image_multiply_gaussian( uint3, float, cuFloatComplex*, unsigned int );
template void image_multiply_gaussian( uint4, float, cuFloatComplex*, unsigned int );

template void clear_image( unsigned int, cuFloatComplex val, cuFloatComplex* );
template void clear_image( unsigned int, float val, float* );

template void square_image( unsigned int, cuFloatComplex* );
template void square_image( unsigned int, float* );

template void add_images( unsigned int, cuFloatComplex*, cuFloatComplex*, cuFloatComplex* );
template void add_images( unsigned int, float*, float*, float* );

template void add_images( unsigned int, cuFloatComplex*, cuFloatComplex*, unsigned int );
template void add_images( unsigned int, float*, float*, unsigned int );

template void subtract_images( unsigned int, cuFloatComplex*, cuFloatComplex*, cuFloatComplex* );
template void subtract_images( unsigned int, float*, float*, float* );

template void multiply_images( unsigned int, cuFloatComplex*, cuFloatComplex*, cuFloatComplex* );
template void multiply_images( unsigned int, float*, float*, float* );

template void image_copy_with_zero_fill( uint2, uint2, cuFloatComplex*, cuFloatComplex*, unsigned int );
template void image_copy_with_zero_fill( uint3, uint3, cuFloatComplex*, cuFloatComplex*, unsigned int );
template void image_copy_with_zero_fill( uint4, uint4, cuFloatComplex*, cuFloatComplex*, unsigned int );

