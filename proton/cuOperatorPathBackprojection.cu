#include "cuOperatorPathBackprojection.h"
#include "vector_td_utilities.h"
#include "vector_td_io.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "check_CUDA.h"

#include <vector>

#include <stdio.h>
#include "hoNDArray_fileio.h"

#include "proton_kernels.cu"


#define MAX_THREADS_PER_BLOCK 128
#define BLOCKS_PER_GRID 65535
#define MAX_BLOCKS 4096*4

/*template <typename T> __inline__ __host__ __device__ int sgn(T val)
{
    return (val > T(0)) - (val < T(0));
}*/


template<class REAL> void cuOperatorPathBackprojection<REAL>
    ::mult_M( cuNDArray<REAL>* in, cuNDArray<REAL>* out, bool accumulate ) {
	 if( !in || !out){
	    BOOST_THROW_EXCEPTION(runtime_error( "cuOperatorPathBackprojection: mult_M empty data pointer"));;
	  }
	 if (!accumulate) clear(out);

	  int dims =  out->get_number_of_elements();

	  int threadsPerBlock =std::min(dims,MAX_THREADS_PER_BLOCK);
		dim3 dimBlock( threadsPerBlock);
		int totalBlocksPerGrid = (dims+threadsPerBlock-1)/threadsPerBlock;
		dim3 dimGrid(std::min(totalBlocksPerGrid,MAX_BLOCKS));
		typename uintd<3>::Type _dims = vector_to_uintd<3>( *(in->get_dimensions().get()) );

		// Invoke kernel
		int batchSize = dimGrid.x*dimBlock.x;
		//std::cout << "Starting forward kernel with grid " << dimGrid.x << " " << dimGrid.y << " " << dimGrid.z << std::endl;
		cudaFuncSetCacheConfig(forward_kernel<REAL>, cudaFuncCachePreferL1);
		 for (int offset = 0; offset < (dims+batchSize); offset += batchSize){
		  forward_kernel<REAL><<< dimGrid, dimBlock >>> (in->get_data_ptr(), out->get_data_ptr(),splines->get_data_ptr(),physical_dims, _dims, dims,offset);
	  }

	  cudaDeviceSynchronize();
	  CHECK_FOR_CUDA_ERROR();
	  if (this->weights.get()){
	  	*out *= *this->weights;
	  }

}

template<class REAL> void cuOperatorPathBackprojection<REAL>
    ::mult_MH( cuNDArray<REAL>* in, cuNDArray<REAL>* out, bool accumulate ) {
	 if( !in || !out){
		 BOOST_THROW_EXCEPTION(runtime_error("cuOperatorPathBackprojection: mult_MH empty data pointer"));
	  }
	 if (!accumulate) clear(out);

	 int dims =  in->get_number_of_elements();
	 int threadsPerBlock =std::min(dims,MAX_THREADS_PER_BLOCK);
	 dim3 dimBlock( threadsPerBlock);
	 int totalBlocksPerGrid = (dims+threadsPerBlock-1)/threadsPerBlock;
	 dim3 dimGrid(std::min(totalBlocksPerGrid,MAX_BLOCKS));
	 typename uintd<3>::Type _dims = vector_to_uintd<3>( *(out->get_dimensions().get()) );

	 // Invoke kernel
	 int batchSize = dimBlock.x*dimGrid.x;
	 //std::cout << "Starting forward kernel with grid " << dimGrid.x << " " << dimGrid.y << " " << dimGrid.z << std::endl;

	 std::cout << "Dimensions: " << _dims << std::endl;
	 std::cout << "Elements: " << out->get_number_of_elements() << " "  << prod(_dims) << std::endl;
	 cudaFuncSetCacheConfig(backwards_kernel<REAL>, cudaFuncCachePreferL1);
	 for (int offset = 0; offset <  (dims+batchSize); offset += batchSize){
		 backwards_kernel<REAL><<< dimGrid, dimBlock >>> (in->get_data_ptr(), out->get_data_ptr(),splines->get_data_ptr(),physical_dims, _dims, dims,offset);
	 }

	 cudaDeviceSynchronize();
	 CHECK_FOR_CUDA_ERROR();
}

template<class REAL> void cuOperatorPathBackprojection<REAL>
::mult_MH_M( cuNDArray<REAL>* in, cuNDArray<REAL>* out, bool accumulate ) {

	cuNDArray<REAL> tmp;

	std::vector<unsigned int> tmp_dim = *(splines->get_dimensions().get());
	tmp_dim[0] /= 4;

	tmp.create(&tmp_dim);

	mult_M(in,&tmp);

	mult_MH(&tmp,out);

}
template<class REAL> void cuOperatorPathBackprojection<REAL>
::setup(boost::shared_ptr< cuNDArray< vector_td<REAL,3> > > splines,  vector_td<REAL,3> physical_dims,  boost::shared_ptr< cuNDArray< REAL > > projections,  boost::shared_ptr< cuNDArray<REAL> > weights,vector_td<REAL,3> origin, REAL background){
	this->weights=weights;
	setup(splines,physical_dims,projections,origin, background);
}
template<class REAL> void cuOperatorPathBackprojection<REAL>
::setup(boost::shared_ptr< cuNDArray< vector_td<REAL,3> > > splines,  vector_td<REAL,3> physical_dims,boost::shared_ptr< cuNDArray< REAL > > projections, vector_td<REAL,3> origin, REAL background){

	this->splines = splines;
	this->physical_dims = physical_dims;
	this->background = background;
	this->origin = origin;

	int dims = splines->get_number_of_elements()/4;


	if (!this->splines->get_data_ptr()) BOOST_THROW_EXCEPTION(runtime_error("Splines data is empty."));
	if (!projections->get_data_ptr()) BOOST_THROW_EXCEPTION(runtime_error("Projections data is empty."));
	if (projections->get_number_of_elements() != dims) BOOST_THROW_EXCEPTION(runtime_error("Projections data does not match splines."));
	 

	 int threadsPerBlock =std::min(dims,MAX_THREADS_PER_BLOCK);
	 dim3 dimBlock( threadsPerBlock);

	 int totalBlocksPerGrid = (dims+threadsPerBlock-1)/threadsPerBlock;
	 dim3 dimGrid(std::min(totalBlocksPerGrid,MAX_BLOCKS));



	 int batchSize = dimGrid.x*dimBlock.x;

	 for (int offset = 0; offset <  (dims+batchSize); offset += batchSize){

	   crop_splines_kernel<<< dimGrid, dimBlock >>> (this->splines->get_data_ptr(),projections->get_data_ptr(),physical_dims,origin,dims,background,offset);
	   cudaThreadSynchronize();
	   CHECK_FOR_CUDA_ERROR();
	 }
	 cudaThreadSynchronize();
   CHECK_FOR_CUDA_ERROR();
	 if (rescale_dirs){
		 for (int offset = 0; offset <  (dims+batchSize); offset += batchSize){
		   rescale_directions_kernel<<< dimGrid, dimBlock >>> (this->splines->get_data_ptr(),projections->get_data_ptr(),physical_dims,(int) dims,offset);
		   CHECK_FOR_CUDA_ERROR();
		 }

	 }
	 cudaThreadSynchronize();

	 for (int offset = 0; offset <  (dims+batchSize); offset += batchSize){
	   points_to_coefficients<<< dimGrid, dimBlock >>>(this->splines->get_data_ptr(),dims,offset);
	   CHECK_FOR_CUDA_ERROR();
	 }



}
// Instantiations
template class cuOperatorPathBackprojection<float>;
