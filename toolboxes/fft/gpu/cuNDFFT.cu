#include "cuNDFFT.h"
#include "cudaDeviceManager.h"

using namespace Gadgetron;

template<class T> __global__ void timeswitch_kernel(T* data, int dimsize, int batchsize, size_t nelements){
	int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < nelements){
		int index = (idx/batchsize)%dimsize;
		if (index & 1) //Check if number is odd
			data[idx] *= -1;
	}
}

template<class T> __global__ void timeswitch_kernel1D(T* data, size_t nelements){

	int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < nelements){
		data[idx] *= (-int(threadIdx.x & 1)*2+1); //Multiply by -1 if x  is odd

	}
}


template<class T> __global__ void timeswitch_kernel2D(T* data, size_t nelements){

	int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < nelements)
		data[idx] *= (-int(threadIdx.x & 1)*2+1)*int(-(blockIdx.x & 1)*2+1); //Multiply by -1 if x or y coordinate is odd, but not if both are
}

template<class T> __global__ void timeswitch_kernel3D(T* data, size_t nelements){
	int idx = ((blockIdx.z*gridDim.y+blockIdx.y)*gridDim.x + blockIdx.x)*blockDim.x + threadIdx.x;
	if (idx < nelements)
		data[idx] *= (-int(threadIdx.x & 1)*2+1)*(-int(blockIdx.x & 1)*2+1)*(-int(blockIdx.y & 1)*2+1); //Multiply by -1 if x or y coordinate is odd, but not if both are
}

template<class T> void Gadgetron::timeswitch1D(cuNDArray<complext<T> >* inout){
	if (inout->get_size(0) > cudaDeviceManager::Instance()->max_blockdim())
		return timeswitch(inout,0);
	dim3 dimBlock(inout->get_size(0));
	size_t max_grid = cudaDeviceManager::Instance()->max_griddim();
	size_t nelements = inout->get_number_of_elements();
	size_t gridX = std::max(std::min(nelements/dimBlock.x,max_grid),size_t(1));
	size_t gridY = std::max(size_t(1),nelements/(gridX*dimBlock.x));
	dim3 dimGrid(gridX,gridY);
	timeswitch_kernel1D<<<dimGrid,dimBlock>>>(inout->get_data_ptr(),nelements);
}

template<class T> void Gadgetron::timeswitch2D(cuNDArray<complext<T> >* inout){
	if (inout->get_size(0) > cudaDeviceManager::Instance()->max_blockdim()){
		timeswitch(inout,0);
		timeswitch(inout,1);
		return;
	}
	dim3 dimBlock(inout->get_size(0));
	size_t max_grid = cudaDeviceManager::Instance()->max_griddim();
	size_t nelements = inout->get_number_of_elements();
	size_t gridX = inout->get_size(1);
	size_t gridY = std::max(size_t(1),nelements/(gridX*dimBlock.x));
	dim3 dimGrid(gridX,gridY);
	timeswitch_kernel2D<<<dimGrid,dimBlock>>>(inout->get_data_ptr(),nelements);

}


template<class T> void Gadgetron::timeswitch3D(cuNDArray<complext<T> >* inout){
	if (inout->get_size(0) > cudaDeviceManager::Instance()->max_blockdim()){
		timeswitch(inout,0);
		timeswitch(inout,1);
		timeswitch(inout,2);
		return;
	}
	dim3 dimBlock(inout->get_size(0));
	size_t max_grid = cudaDeviceManager::Instance()->max_griddim();
	size_t nelements = inout->get_number_of_elements();
	size_t gridX = inout->get_size(1);
	size_t gridY = inout->get_size(2);
	size_t gridZ = std::max(size_t(1),nelements/(gridX*dimBlock.x*gridY));
	dim3 dimGrid(gridX,gridY,gridZ);
	timeswitch_kernel3D<<<dimGrid,dimBlock>>>(inout->get_data_ptr(),nelements);
}

template<class T> void Gadgetron::timeswitch(cuNDArray<complext<T> >* inout, int dim_to_transform){


	size_t batchsize = 1;
	for (int i = 0; i < dim_to_transform; i++)
		batchsize *= inout->get_size(i);

	size_t dimsize = inout->get_size(dim_to_transform);

	size_t nelements = inout->get_number_of_elements();

	size_t max_block = cudaDeviceManager::Instance()->max_blockdim();
	dim3 dimBlock(std::min(max_block,nelements));

	size_t max_grid = cudaDeviceManager::Instance()->max_griddim();
	size_t gridX = std::max(std::min(nelements/dimBlock.x,max_grid),size_t(1));
	size_t gridY = std::max(size_t(1),nelements/(gridX*dimBlock.x));

	dim3 dimGrid(gridX,gridY);

	timeswitch_kernel<<<dimGrid,dimBlock>>>(inout->get_data_ptr(),dimsize,batchsize,nelements);



}



template void Gadgetron::timeswitch<float>(cuNDArray<float_complext>*, int);
template void Gadgetron::timeswitch<double>(cuNDArray<double_complext>*, int);

template void Gadgetron::timeswitch1D<float>(cuNDArray<float_complext>*);
template void Gadgetron::timeswitch1D<double>(cuNDArray<double_complext>*);

template void Gadgetron::timeswitch2D<float>(cuNDArray<float_complext>*);
template void Gadgetron::timeswitch2D<double>(cuNDArray<double_complext>*);

template void Gadgetron::timeswitch3D<float>(cuNDArray<float_complext>*);
template void Gadgetron::timeswitch3D<double>(cuNDArray<double_complext>*);
