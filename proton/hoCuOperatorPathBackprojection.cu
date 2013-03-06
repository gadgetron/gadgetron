#include "hoCuOperatorPathBackprojection.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

#include "check_CUDA.h"

#include <vector>

#include <stdio.h>
#include <stdexcept>
#include <algorithm>



#include "proton_kernels.cu"


#define MAX_THREADS_PER_BLOCK 128
#define BLOCKS_PER_GRID 65535
#define MAX_BLOCKS 4096*4


template<class REAL> size_t hoCuOperatorPathBackprojection<REAL>::calculate_batch_size(){

	int mem_per_proton = 13*sizeof(REAL); // Need 12 REALS for the splines and 1 for the projection
	size_t free;
	size_t total;

	int res = cudaMemGetInfo(&free,&total);

	//printReturn(res);

	std::cout << "Free memory in MB: " << free << std::endl;
	return 1024*1024*(free/(1024*1024*mem_per_proton)); //Divisons by 1024*1024 to ensure MB batch size
}

template<class REAL> void hoCuOperatorPathBackprojection<REAL>
    ::mult_M( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate ) {
	 if( !in || !out){
	    BOOST_THROW_EXCEPTION(runtime_error("cuOperatorPathBackprojection: mult_M empty data pointer"));
	  }

	 cuNDArray<REAL> image(in);
	 size_t max_batch_size = calculate_batch_size();
	 size_t elements = out->get_number_of_elements();
	 size_t offset = 0;


	 for (size_t n = 0; n < (elements+max_batch_size-1)/max_batch_size; n++){

		 size_t batch_size = std::min(max_batch_size,elements-offset);
		 std::vector<unsigned int> batch_dim;
		 batch_dim.push_back(batch_size);

		 hoCuNDArray<REAL> out_view(&batch_dim,out->get_data_ptr()+offset); //This creates a "view" of out

		 cuNDArray<REAL> out_dev(&out_view);
		 batch_dim.push_back(4);
		 hoCuNDArray<vector_td<REAL,3> > splines_view(&batch_dim,splines->get_data_ptr()+offset*4); // This creates a "view" of splines
		 cuNDArray<vector_td<REAL,3> > splines_dev(&splines_view);
		 if (!accumulate) cudaMemset(out_dev.get_data_ptr(),0,sizeof(REAL)*batch_size);

		 int threadsPerBlock = std::min((int)batch_size,MAX_THREADS_PER_BLOCK);
		 dim3 dimBlock( threadsPerBlock);
		 int totalBlocksPerGrid = (batch_size+threadsPerBlock-1)/threadsPerBlock;
		 dim3 dimGrid(std::min(totalBlocksPerGrid,MAX_BLOCKS));
		 typename uintd<3>::Type _dims = vector_to_uintd<3>( *(image.get_dimensions().get()) );

		 // Invoke kernel
		 int offset_k = 0;
		 //std::cout << "Starting forward kernel with grid " << dimGrid.x << " " << dimGrid.y << " " << dimGrid.z << std::endl;
		 for (int i = 0; i <= (totalBlocksPerGrid+dimGrid.x-1)/dimGrid.x; i++){
			forward_kernel<REAL><<< dimGrid, dimBlock >>> (image.get_data_ptr(), out_dev.get_data_ptr(),splines_dev.get_data_ptr(),physical_dims, _dims, batch_size,offset_k);
			offset_k += dimBlock.x*dimGrid.x;
		 }
		 //cudaDeviceSynchronize();
		 CHECK_FOR_CUDA_ERROR();
		 if (this->weights.get()){
			//cuNDA_scale(this->weights.get(),&out_dev);
			throw std::runtime_error("Weights currently not supported");
		 }

		 cudaMemcpy(out_view.get_data_ptr(),out_dev.get_data_ptr(),batch_size*sizeof(REAL),cudaMemcpyDeviceToHost); //Copies back the data to the host

		offset += batch_size;

	 }
}

template<class REAL> void hoCuOperatorPathBackprojection<REAL>
    ::mult_MH( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate ) {
	 if( !in || !out){
	    BOOST_THROW_EXCEPTION(runtime_error("cuOperatorPathBackprojection: mult_MH empty data pointer"));
	  }
	 if (this->weights.get()){
		 BOOST_THROW_EXCEPTION(runtime_error("Weights unsupported at the moment!"));
	 }

	 cuNDArray<REAL> image(out);
	 CHECK_FOR_CUDA_ERROR();
	 if (!accumulate) cudaMemset(image.get_data_ptr(),0,image.get_number_of_elements()*sizeof(REAL));
	 CHECK_FOR_CUDA_ERROR();
	 size_t max_batch_size = calculate_batch_size();
	 size_t elements = in->get_number_of_elements();
	 size_t offset = 0;

	 for (size_t n = 0; n < (elements+max_batch_size-1)/max_batch_size; n++){
		 size_t batch_size = std::min(max_batch_size,elements-offset);
		 std::vector<unsigned int> batch_dim;
		 batch_dim.push_back(batch_size);

		 hoCuNDArray<REAL> in_view(&batch_dim,in->get_data_ptr()+offset); //This creates a "view" of in

		 cuNDArray<REAL> in_dev(&in_view);
		 batch_dim.push_back(4);
		 hoCuNDArray<vector_td<REAL,3> > splines_view(&batch_dim,splines->get_data_ptr()+offset*4); // This creates a "view" of splines
		 cuNDArray<vector_td<REAL,3> > splines_dev(&splines_view);


		 int threadsPerBlock = std::min((int)batch_size,MAX_THREADS_PER_BLOCK);
		 dim3 dimBlock( threadsPerBlock);
		 int totalBlocksPerGrid = (batch_size+threadsPerBlock-1)/threadsPerBlock;
		 dim3 dimGrid(std::min(totalBlocksPerGrid,MAX_BLOCKS));
		typename uintd<3>::Type _dims = vector_to_uintd<3>( *(image.get_dimensions().get()) );

		// Invoke kernel
		int offset_k = 0;
		//std::cout << "Starting forward kernel with grid " << dimGrid.x << " " << dimGrid.y << " " << dimGrid.z << std::endl;
		for (int i = 0; i <= (totalBlocksPerGrid+dimGrid.x-1)/dimGrid.x; i++){
			backwards_kernel<REAL><<< dimGrid, dimBlock >>> (in_dev.get_data_ptr(), image.get_data_ptr(),splines_dev.get_data_ptr(),physical_dims, _dims, batch_size, offset_k);
			offset_k += dimBlock.x*dimGrid.x;
		}


			CHECK_FOR_CUDA_ERROR();
			offset += batch_size;
	 }

	 cudaMemcpy(out->get_data_ptr(),image.get_data_ptr(),image.get_number_of_elements()*sizeof(REAL),cudaMemcpyDeviceToHost);
}

template<class REAL> void hoCuOperatorPathBackprojection<REAL>
::mult_MH_M( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate ) {

	hoCuNDArray<REAL> tmp;

	std::vector<unsigned int> tmp_dim = *(splines->get_dimensions().get());
	tmp_dim[0] /= 4;

	tmp.create(&tmp_dim);

	mult_M(in,&tmp);

	mult_MH(&tmp,out);

}
template<class REAL> void hoCuOperatorPathBackprojection<REAL>
::setup(boost::shared_ptr< hoCuNDArray< vector_td<REAL,3> > > splines,  vector_td<REAL,3> physical_dims,  boost::shared_ptr< hoCuNDArray< REAL > > projections,  boost::shared_ptr< hoCuNDArray<REAL> > weights, vector_td<REAL,3> origin,REAL background){
	this->weights=weights;
	setup(splines,physical_dims,projections,origin,background);
}
template<class REAL> void hoCuOperatorPathBackprojection<REAL>
::setup(boost::shared_ptr< hoCuNDArray< vector_td<REAL,3> > > splines,  vector_td<REAL,3> physical_dims,boost::shared_ptr< hoCuNDArray< REAL > > projections, vector_td<REAL,3> origin, REAL background){

	this->splines = splines;
	this->physical_dims = physical_dims;
	this->background = background;
	this->origin = origin;
	std::cout  << "Calling setup function nowish" << std::endl;
	std::vector<unsigned int> dims = *(splines->get_dimensions().get());
	dims[0] /= 4;
	size_t max_batch_size = calculate_batch_size();


	size_t elements = splines->get_number_of_elements()/4;
	size_t offset = 0;
	for (size_t n = 0; n < (elements-1)/max_batch_size+1; n++){


		size_t batch_size = std::min(max_batch_size,elements-offset);
		std::vector<unsigned int> batch_dim;
		batch_dim.push_back(batch_size);

		hoCuNDArray<REAL> projections_view(&batch_dim,projections->get_data_ptr()+offset); //This creates a "view" of out


		cuNDArray<REAL> projections_dev(&projections_view);
		CHECK_FOR_CUDA_ERROR();
		batch_dim.push_back(4);
		hoCuNDArray<vector_td<REAL,3> > splines_view(&batch_dim,splines->get_data_ptr()+offset*4); // This creates a "view" of splines
		cuNDArray<vector_td<REAL,3> > splines_dev(&splines_view);
		 CHECK_FOR_CUDA_ERROR();




		int threadsPerBlock = std::min((int)batch_size, MAX_THREADS_PER_BLOCK);
		dim3 dimBlock( threadsPerBlock);
		int totalBlocksPerGrid = (batch_size-1)/MAX_THREADS_PER_BLOCK+1;
		dim3 dimGrid(std::min(totalBlocksPerGrid,MAX_BLOCKS));


		// Invoke kernel
		int offset_k = 0;
		//std::cout << "Starting forward kernel with grid " << dimGrid.x << " " << dimGrid.y << " " << dimGrid.z << std::endl;
		for (int i = 0; i <= (totalBlocksPerGrid-1)/MAX_BLOCKS+1; i++){
			 crop_splines_kernel<<< dimGrid, dimBlock >>> (splines_dev.get_data_ptr(),projections_dev.get_data_ptr(),physical_dims, origin, batch_size,background,offset_k);
			 offset_k += dimGrid.x*dimBlock.x;

		}
		cudaThreadSynchronize();
		// Invoke kernel
		if (rescale_dirs){
			offset_k = 0;
			//std::cout << "Starting forward kernel with grid " << dimGrid.x << " " << dimGrid.y << " " << dimGrid.z << std::endl;
			for (int i = 0; i <= (totalBlocksPerGrid-1)/MAX_BLOCKS+1; i++){
				 rescale_directions_kernel<<< dimGrid, dimBlock >>> (splines_dev.get_data_ptr(),projections_dev.get_data_ptr(),physical_dims, batch_size,offset_k);
				 offset_k += dimGrid.x*dimBlock.x;

			}
			cudaThreadSynchronize();
		}
	 offset_k = 0;
	 for (int i = 0; i <= (totalBlocksPerGrid-1)/MAX_BLOCKS+1; i++){
	   points_to_coefficients<<< dimGrid, dimBlock >>>(splines_dev.get_data_ptr(),batch_size,offset_k);
	   offset_k += dimGrid.x*dimBlock.x;
	   CHECK_FOR_CUDA_ERROR();
	 }

	 cudaMemcpy(splines_view.get_data_ptr(),splines_dev.get_data_ptr(),batch_size*4*sizeof(vector_td<REAL,3>),cudaMemcpyDeviceToHost);
	 cudaMemcpy(projections_view.get_data_ptr(),projections_dev.get_data_ptr(),batch_size*sizeof(REAL),cudaMemcpyDeviceToHost);
	 offset += batch_size;

	 }


}
// Instantiations
template class hoCuOperatorPathBackprojection<float>;
