#include <cublas.h>
#include <cufft.h>

#include "batchfft.hcu"
#include "uint_util.hcu"

/**
Kernel for performing the FFT wrap

*/
template<class T> __global__ void fft_shift_kernel( cuFloatComplex* data_in, cuFloatComplex* data_out, T dim, unsigned int num_images)
{

	// This is the current pixel number
	unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	// Number of elements
	unsigned int num_elements = prod(dim);

	// Do nothing if index is out of range
	if( idx < num_elements ){

		//Lets get the coordinates for this pixel
		T src_co = idx_to_co(idx, dim);

		//Where should this data go?
		T dst_co = (src_co+(dim>>1))%dim;

		for( unsigned int i=0; i<num_images; i++ ){
			unsigned int offset = i*num_elements;
			data_out[co_to_idx(dst_co,dim)+offset] = data_in[idx+offset];
		}  
	}
}

/**
Kernel for performing the FFT wrap

This one also permutes the dim_to_traf dimension to the first dimensions for a 1D FFT

*/
template<class T> __global__ void fft_shift_permute_kernel( unsigned int num_elements, cuFloatComplex* data_in, cuFloatComplex* data_out, T dim, T offset, unsigned int dim_to_traf, unsigned int d )
{	
	//This is the current pixel number
	unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

	// Do nothing if index is out of range
	if( idx < num_elements ){

		//Lets get the coordinates for this pixel
		T src_co = idx_to_co(idx, dim)+offset;

		//Where should this data go?
		T dst_co = (src_co+(dim>>1))%dim;
		T new_dim;

		// 'd' is the "direction": 0 or 1, i.e. shift before or after the FFT
		if (d == 1)
		{
			dst_co = shift_down(dst_co,dim_to_traf);
			new_dim = shift_down(dim,dim_to_traf);
		} else {
			dst_co = shift_up(dst_co,dim_to_traf);
			new_dim = shift_up(dim,dim_to_traf);
		}

		data_out[co_to_idx(dst_co,new_dim)] = data_in[idx];
	}
}

template<class T> void 
fft_shift(cuFloatComplex* data_in, cuFloatComplex* data_out, T dim, unsigned int num_images )
{  
	if( data_in == data_out ){
		printf("\nError in fft shifter. Input and output pointer cannot overlap.\n");
		return;
	}

	dim3 blockDim(512,1,1);
	dim3 gridDim((int)ceil((double)prod(dim)/(double)blockDim.x),1,1);

	fft_shift_kernel<T><<< gridDim, blockDim >>>(data_in, data_out, dim, num_images);

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		printf("\nCuda error detected in 'fft_shift': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
		exit(1);
	}
}


/**
This function will perform the FFT along the dimension specified by dim_to_trafo
in the direction (CUFFT_FORWARD or CUFFT_INVERSE) specified by direction

*/
template<class T> __host__ bool ft_1d(cuFloatComplex* data, T dim, T offset, unsigned int dim_to_traf, int direction, bool do_scale, bool do_shift )
{
	if( dim_to_traf != 1 && !do_shift )
		printf("\nWARNING: 1D FFT of multidimensional dataset WITH NO SHIFT is requested. This is a problem if data is not already permuted!\n");

	unsigned int Nx = ((unsigned int*)&dim)[dim_to_traf];

	cuFloatComplex* temp;  
	dim3 blockDim(512,1, 1);
	dim3 gridDim((int)ceil((double)prod(dim)/(double)blockDim.x),1,1);

	if( do_shift ){
		cudaMalloc( (void **) &temp,sizeof(cuFloatComplex)*prod(dim));
		fft_shift_permute_kernel<T><<< gridDim, blockDim >>>(prod(dim), data, (cuFloatComplex*)temp, dim, offset, dim_to_traf, 1 );

		cudaError_t err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected in 'fft_shift_permute_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}
	}  
	else
		temp = (cuFloatComplex*) data;

	cufftHandle plan;
	cufftResult res;
	
	res = cufftPlan1d(&plan, Nx, CUFFT_C2C, prod(dim)/Nx); // Batched 1D FFTs

	if( res != CUFFT_SUCCESS ){
		printf("\nFATAL ERROR in 'cufftPlan2d': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
		exit(1);
	}

	res = cufftExecC2C(plan, temp, temp, direction);

	if( res != CUFFT_SUCCESS ){
		printf("\nFATAL ERROR in 'cufftExecC2C': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
		exit(1);
	}

	if (direction == CUFFT_INVERSE && do_scale)
	{
		cublasCscal (prod(dim), make_cuFloatComplex(1.0f/Nx,0.0f), (cuFloatComplex*)temp, 1);
	}

	if( do_shift ){
		fft_shift_permute_kernel<T><<< gridDim, blockDim >>>(prod(dim), (cuFloatComplex*)temp, data, shift_down(dim,dim_to_traf), offset, dim_to_traf, 0 );
	
		cudaError_t err = cudaGetLastError();
		if( err != cudaSuccess ){
			printf("\nCuda error detected in 'fft_shift_permute_kernel': %s. Quitting.\n", cudaGetErrorString(err) ); fflush(stdout);
			exit(1);
		}

		cudaFree(temp);
	}

	res = cufftDestroy(plan);

	if( res != CUFFT_SUCCESS ){
		printf("\nFATAL ERROR in 'cufftDestroy': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
		exit(1);
	}

	return true;
}

/**
This function determines the correct conversion of the dimension parameters based on the dimensionaly of the problem,
it then calls the ft_1d_wrapper.
*/
bool ft_1d_wrapper(cuFloatComplex* data, uint4 dim, unsigned int dim_to_trans, int direction, bool do_scale, bool do_shift )
{
	uint4 offset = make_uint4(0,0,0,0);

	if (dim.z == 1 && dim.w == 1)
	{
		//One or Two dimensionsal dataset
		ft_1d(data, uint4_to_uint2(dim), uint4_to_uint2(offset), dim_to_trans, direction, do_scale, do_shift );
	}
	else if (dim.w == 1)
	{
		//Three dimensional dataset
		ft_1d(data, uint4_to_uint3(dim),uint4_to_uint3(offset), dim_to_trans, direction, do_scale, do_shift );
	}
	else
	{
		//Four dimensional dataset
		ft_1d(data, dim, offset, dim_to_trans, direction, do_scale, do_shift );
	}

	return true;
}

/**
This wrapper is used for a 2D FFT where both dimensions are transformed.

*/
bool ft_2d_wrapper(cuFloatComplex* data, uint2 dim, int direction, unsigned int num_images, bool do_scale, bool shift)
{
	cuFloatComplex* temp;

	if( shift ){
		cudaMalloc( (void **) &temp, sizeof(cuFloatComplex)*prod(dim)*num_images );
		fft_shift<uint2>( data, temp, dim, num_images );
	}
	else
		temp = (cuFloatComplex*)data;

	if( num_images == 1 ){

		// Single FFT

		cufftHandle plan;
		cufftResult res;
		
		res = cufftPlan2d(&plan, dim.y, dim.x, CUFFT_C2C);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'cufftPlan2d': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}

		res = cufftExecC2C(plan, temp, temp, direction);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'cufftExecC2C': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}

		res = cufftDestroy(plan);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'cufftDestroy': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}
	}
	else{

		// Batched 2D FFTs

		batchfftHandle batchplan;
		cufftResult res;
		
		res = batchfftPlan2d(&batchplan, dim.y, dim.x, CUFFT_C2C, num_images);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'batchfftPlan2d': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}

		res = batchfftExecute(batchplan, temp, temp, direction);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'batchfftExecute': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}

		res = batchfftDestroy(&batchplan);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'batchfftDestroy': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}
	}

	if (direction == CUFFT_INVERSE && do_scale)
	{
		cublasCscal( prod(dim)*num_images, make_cuFloatComplex(1.0f/prod(dim),0.0f), (cuFloatComplex*)temp, 1 );
	}

	if( shift ){
		fft_shift<uint2>( (cuFloatComplex*)temp, data, dim, num_images );
		cudaFree( temp );
	}

	return true;
}

/** 
This function is used for a 3D FFT where all dimensions are transformed

*/
bool ft_3d_wrapper(cuFloatComplex* data, uint3 dim, int direction, unsigned int num_images, bool do_scale, bool shift)
{
	cuFloatComplex* temp;

	if(shift){
		cudaMalloc( (void **) &temp, sizeof(cuFloatComplex)*prod(dim)*num_images );
		fft_shift<uint3>(data, (cuFloatComplex*)temp, dim, num_images );
	}
	else
		temp = (cuFloatComplex*)data;

	if( num_images == 1 ){

		// Single FFT

		cufftHandle plan;
		cufftResult res;
		
		res = cufftPlan3d(&plan, dim.z, dim.y, dim.x, CUFFT_C2C);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'cufftPlan3d': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}

		res = cufftExecC2C(plan, temp, temp, direction);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'cufftExecC2C': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}

		res = cufftDestroy(plan);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'cufftDestroy': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}
	}
	else{

		// Batched 3D FFTs
		batchfftHandle batchplan;
		cufftResult res;
		
		res = batchfftPlan3d(&batchplan, dim.z, dim.y, dim.x, CUFFT_C2C, num_images);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'batchfftPlan3d': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}

		res = batchfftExecute(batchplan, temp, temp, direction);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'batchfftExecute': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}

		res = batchfftDestroy(&batchplan);

		if( res != CUFFT_SUCCESS ){
			printf("\nFATAL ERROR in 'batchfftDestroy': %s. Quitting.\n", cudaGetErrorString(cudaGetLastError()));
			exit(1);
		}
	}

	if (direction == CUFFT_INVERSE && do_scale)
	{
		cublasCscal (prod(dim)*num_images, make_cuFloatComplex(1.0f/prod(dim),0.0f), (cuFloatComplex*)temp, 1);
	}

	if(shift){
		fft_shift<uint3>( (cuFloatComplex*)temp, data, dim, num_images );
		cudaFree( temp );
	}


	return true;
}

/**
This is the exported function, which can be called from the main program.
It calls the wrapper.

It performs the FFT from image space to k-space (i.e. forward transform)

*/
__host__ bool I2K(cuFloatComplex* data, uint4 dim, unsigned int dim_to_trans, bool do_scale, bool do_shift )
{
	return ft_1d_wrapper(data,dim,dim_to_trans, CUFFT_FORWARD, do_scale, do_shift );
}


/**
This is the exported function, which can be called from the main program.
It calls the wrapper.

It performs the INVERSE FFT from k-space to image space (inverse transform)

*/
__host__ bool K2I(cuFloatComplex* data, uint4 dim, unsigned int dim_to_trans, bool do_scale, bool do_shift )
{
	return ft_1d_wrapper(data,dim,dim_to_trans, CUFFT_INVERSE, do_scale, do_shift );
}

/**
Transform of all relevant dimensions from image to k space.

Calls the appropriate wrapper.

*/
__host__ bool I2K_ALL(cuFloatComplex* data, uint4 dim, unsigned int num_images, bool do_scale, bool do_shift)
{
	if (dim.z == 1 && dim.w == 1)
	{
		//One or Two dimensionsal dataset
		ft_2d_wrapper(data, uint4_to_uint2(dim), CUFFT_FORWARD, num_images, do_scale, do_shift);
	}
	else if (dim.w == 1)
	{
		//Three dimensional dataset
		ft_3d_wrapper(data, uint4_to_uint3(dim), CUFFT_FORWARD, num_images, do_scale, do_shift);
	}
	else
	{
		if( num_images>1 ){
			printf("\nI2K_ALL cannot handles batches at the moment!. Quitting.\n");
			exit(1);
		}

		//Four dimensional dataset
		//The CUFFT library doesn't support 4D transforms, so we'll have to hack it
		for (unsigned int i = 0; i < dim.w; i++)
		{
			ft_3d_wrapper( (data+i*prod(uint4_to_uint3(dim))), uint4_to_uint3(dim), CUFFT_FORWARD, num_images, do_scale, do_shift );
		}
		ft_1d_wrapper( data,dim,3,CUFFT_FORWARD, do_scale, do_shift ); //Last dimension manually
	}

	return true;
}

__host__ bool I2K_ALL(cuFloatComplex* data, uint3 dim, unsigned int num_images, bool do_scale, bool do_shift)
{
	if (dim.z == 1 )
	{
		//One or Two dimensionsal dataset
		ft_2d_wrapper(data, uint3_to_uint2(dim),CUFFT_FORWARD, num_images, do_scale, do_shift);
	}
	else
	{
		//Three dimensional dataset
		ft_3d_wrapper(data, dim, CUFFT_FORWARD, num_images, do_scale, do_shift);
	}

	return true;
}

__host__ bool I2K_ALL(cuFloatComplex* data, uint2 dim, unsigned int num_images, bool do_scale, bool do_shift)
{
	//One or Two dimensionsal dataset
	ft_2d_wrapper(data, dim, CUFFT_FORWARD, num_images, do_scale, do_shift);

	return true;
}

/**
Transform of all relevant dimensions from k-space to image space.

Calls the appropriate wrapper.

*/
__host__ bool K2I_ALL(cuFloatComplex* data, uint4 dim, unsigned int num_images, bool do_scale, bool do_shift)
{
	if (dim.z == 1 && dim.w == 1)
	{
		//One or Two dimensionsal dataset
		ft_2d_wrapper(data, uint4_to_uint2(dim),CUFFT_INVERSE, num_images, do_scale, do_shift);
	}
	else if (dim.w == 1)
	{
		//Three dimensional dataset
		ft_3d_wrapper(data, uint4_to_uint3(dim), CUFFT_INVERSE, num_images, do_scale, do_shift);
	}
	else
	{

		if( num_images>1 ){
			printf("\nK2I_ALL cannot handles batches at the moment!. Quitting.\n");
			exit(1);
		}

		//Four dimensional dataset
		//The CUFFT library doesn't support 4D transforms, so we'll have to hack it
		for (unsigned int i = 0; i < dim.w; i++)
		{
			ft_3d_wrapper((data+i*prod(uint4_to_uint3(dim))), uint4_to_uint3(dim), CUFFT_INVERSE, num_images, do_scale, do_shift);
		}
		ft_1d_wrapper( data, dim, 3, CUFFT_INVERSE, do_scale, do_shift ); //Last dimension manually
	}

	return true;
}

__host__ bool K2I_ALL(cuFloatComplex* data, uint3 dim, unsigned int num_images, bool do_scale, bool do_shift)
{
	if (dim.z == 1 )
	{
		//One or Two dimensionsal dataset
		ft_2d_wrapper(data, uint3_to_uint2(dim), CUFFT_INVERSE, num_images, do_scale, do_shift);
	}
	else
	{
		//Three dimensional dataset
		ft_3d_wrapper(data, dim, CUFFT_INVERSE, num_images, do_scale, do_shift);
	}

	return true;
}

__host__ bool K2I_ALL(cuFloatComplex* data, uint2 dim, unsigned int num_images, bool do_scale, bool do_shift)
{
	//One or Two dimensionsal dataset
	ft_2d_wrapper(data, dim, CUFFT_INVERSE, num_images, do_scale, do_shift);

	return true;
}

// template instantiation

template void fft_shift(cuFloatComplex*, cuFloatComplex*, uint2, unsigned int );
template void fft_shift(cuFloatComplex*, cuFloatComplex*, uint3, unsigned int );
template void fft_shift(cuFloatComplex*, cuFloatComplex*, uint4, unsigned int );
