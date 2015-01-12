#include "cuNDFFT.h"
#include "cudaDeviceManager.h"

using namespace Gadgetron;

template<class T> __global__ void timeswitch_kernel(T* data, unsigned int dimsize, unsigned int batchsize, size_t nelements){
	unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
	if (idx < nelements){
		size_t index = (idx/batchsize)%dimsize;
		if (index%2 == 1)
			data[idx] *= -1;
	}
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



template EXPORTGPUFFT void Gadgetron::timeswitch<float>(cuNDArray<float_complext>*, int);

template EXPORTGPUFFT void Gadgetron::timeswitch<double>(cuNDArray<double_complext>*, int);
