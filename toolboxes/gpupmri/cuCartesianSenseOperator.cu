#include "cuCartesianSenseOperator.h"
#include "cuNDFFT.h"
#include "ndarray_vector_td_utilities.h"
#include <sstream>

template<class REAL> __global__ void 
sample_array_kernel( complext<REAL> *in, complext<REAL> *out,
		     unsigned int *idx, 
		     unsigned long image_elements,
		     unsigned long int samples,
		     unsigned int coils )
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < samples) {
    for (unsigned int i = 0; i < coils; i++) {
      out[idx_in + i*samples].vec[0] += in[idx[idx_in] + i*image_elements].vec[0];
      out[idx_in + i*samples].vec[1] += in[idx[idx_in] + i*image_elements].vec[1];
    }
  }
}

template<class REAL> __global__ void 
insert_samples_kernel( complext<REAL> *in, complext<REAL> *out,
				       unsigned int *idx, 
				       unsigned long image_elements,
				       unsigned long int samples,
				       unsigned int coils )
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < samples) {
    for (unsigned int i = 0; i < coils; i++) {
      out[idx[idx_in] + i*image_elements].vec[0] += in[idx_in + i*samples].vec[0];
      out[idx[idx_in] + i*image_elements].vec[1] += in[idx_in + i*samples].vec[1];
    }
  }
}

template<class REAL, unsigned int D> void
cuCartesianSenseOperator<REAL,D>::mult_M( cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate )
{

  int ret = this->_set_device();
  if( ret<0 ){
    throw cuda_error("cuCartesianSenseOperator::mult_M: unable to set device");

  }
  
  if (!(in->dimensions_equal(this->get_domain_dimensions().get())) || !(out->dimensions_equal(this->get_codomain_dimensions().get())) ) {

    throw std::runtime_error("cuCartesianSenseOperator::mult_M dimensions mismatch");

  }


  std::vector<unsigned int> full_dimensions = *this->get_domain_dimensions();
  full_dimensions.push_back(this->ncoils_);
  cuNDArray<_complext> tmp(&full_dimensions);



  mult_csm(in,&tmp);

  cuNDFFT<_complext> ft;
  std::vector<unsigned int> ft_dims;
  for (unsigned int i = 0; i < this->get_domain_dimensions()->size(); i++) {
    ft_dims.push_back(i);
  }

  ft.fft(&tmp, &ft_dims);

  if (!accumulate) 
  	out->clear();


  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)idx_->get_number_of_elements()/blockDim.x), 1, 1 );
  sample_array_kernel<REAL><<< gridDim, blockDim >>>( tmp.get_data_ptr(), out->get_data_ptr(), idx_->get_data_ptr(),
						      in->get_number_of_elements(), idx_->get_number_of_elements(), this->ncoils_);
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
     std::stringstream ss;
     ss <<"cuCartesianSenseOperator::mult_M : Unable to sample data: " <<
      cudaGetErrorString(err);
     throw cuda_error(ss.str());

  }

  this->_restore_device();
  
}

template<class REAL, unsigned int D> void
cuCartesianSenseOperator<REAL,D>::mult_MH(cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate)
{
  this->_set_device();


  if (!(out->dimensions_equal(this->get_domain_dimensions().get())) || 
      !(in->dimensions_equal(this->get_codomain_dimensions().get())) ) {
    throw std::runtime_error( "cuCartesianSenseOperator::mult_MH dimensions mismatch");

  }

  std::vector<unsigned int> tmp_dimensions = *this->get_domain_dimensions();
  tmp_dimensions.push_back(this->ncoils_);

  cuNDArray<_complext> tmp(&tmp_dimensions);

  tmp.clear();

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)idx_->get_number_of_elements()/blockDim.x), 1, 1 );
  insert_samples_kernel<REAL><<< gridDim, blockDim >>>( in->get_data_ptr(), tmp.get_data_ptr(),
							idx_->get_data_ptr(),out->get_number_of_elements(),
							idx_->get_number_of_elements(), this->ncoils_);
  
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
  	std::stringstream ss;
    ss << "cuCartesianSenseOperator::mult_EM : Unable to insert samples into array: " <<
      cudaGetErrorString(err);
    throw cuda_error(ss.str());
  }

  cuNDFFT<_complext> ft;
  std::vector<unsigned int> ft_dims;
  for (unsigned int i = 0; i < this->get_domain_dimensions()->size(); i++) {
    ft_dims.push_back(i);
  }

  ft.ifft(&tmp, &ft_dims);

  if (!accumulate) 
    out->clear();
  
  mult_csm_conj_sum(&tmp,out);

  this->_restore_device();

}

//
// Instantiations
//

template class EXPORTGPUPMRI cuCartesianSenseOperator<float,2>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<float,3>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<float,4>;

template class EXPORTGPUPMRI cuCartesianSenseOperator<double,2>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<double,3>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<double,4>;

