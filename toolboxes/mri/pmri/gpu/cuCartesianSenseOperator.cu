#include "cuCartesianSenseOperator.h"
#include "cuNDFFT.h"

#include <sstream>

using namespace Gadgetron;

template<class REAL> __global__ void 
sample_array_kernel( const complext<REAL> * __restrict__ in, complext<REAL> * __restrict__ out,
		     unsigned int *idx, 
		     unsigned int image_elements,
		     unsigned int samples,
		     unsigned int coils )
{
  unsigned int idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < samples) {
    for (unsigned int i = 0; i < coils; i++) {
      out[idx_in + i*samples]._real += in[idx[idx_in] + i*image_elements]._real;
      out[idx_in + i*samples]._imag += in[idx[idx_in] + i*image_elements]._imag;
    }
  }
}

template<class REAL> __global__ void 
insert_samples_kernel( const complext<REAL> * __restrict__ in, complext<REAL> * __restrict__ out,
		       unsigned int *idx, 
		       unsigned int image_elements,
		       unsigned int samples,
		       unsigned int coils )
{
  unsigned int idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < samples) {
    for (unsigned int i = 0; i < coils; i++) {
      out[idx[idx_in] + i*image_elements]._real += in[idx_in + i*samples]._real;
      out[idx[idx_in] + i*image_elements]._imag += in[idx_in + i*samples]._imag;
    }
  }
}

template<class REAL, unsigned int D> void
cuCartesianSenseOperator<REAL,D>::mult_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate )
{
  if (!(in->dimensions_equal(*this->get_domain_dimensions().get())) || !(out->dimensions_equal(*this->get_codomain_dimensions().get())) ) {
    throw std::runtime_error("cuCartesianSenseOperator::mult_M dimensions mismatch");
  }
  
  std::vector<size_t> full_dimensions = *this->get_domain_dimensions();
  full_dimensions.push_back(this->ncoils_);
  cuNDArray< complext<REAL> > tmp(full_dimensions);

  this->mult_csm(in,&tmp);

  std::vector<size_t> ft_dims;
   for (unsigned int i = 0; i < this->get_domain_dimensions()->size(); i++) {
     ft_dims.push_back(i);
   }

   cuNDFFT<REAL>::instance()->fft(&tmp, &ft_dims);
   
  if (!accumulate) 
    clear(out);

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) std::ceil((double)idx_->get_number_of_elements()/blockDim.x), 1, 1 );
  sample_array_kernel<REAL><<< gridDim, blockDim >>>( tmp.get_data_ptr(), out->get_data_ptr(), idx_->get_data_ptr(),
						      in->get_number_of_elements(), idx_->get_number_of_elements(), this->ncoils_);
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::stringstream ss;
    ss <<"cuCartesianSenseOperator::mult_M : Unable to sample data: " <<
      cudaGetErrorString(err);
    throw cuda_error(ss.str());
  }
}

template<class REAL, unsigned int D> void
cuCartesianSenseOperator<REAL,D>::mult_MH(cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate)
{
  if (!(out->dimensions_equal(this->get_domain_dimensions().get())) || 
      !(in->dimensions_equal(this->get_codomain_dimensions().get())) ) {
    throw std::runtime_error( "cuCartesianSenseOperator::mult_MH dimensions mismatch");

  }

  std::vector<size_t> tmp_dimensions = *this->get_domain_dimensions();
  tmp_dimensions.push_back(this->ncoils_);

  cuNDArray< complext<REAL> > tmp(tmp_dimensions);
  clear(&tmp);

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) std::ceil((double)idx_->get_number_of_elements()/blockDim.x), 1, 1 );
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


  std::vector<size_t> ft_dims;
  for (unsigned int i = 0; i < this->get_domain_dimensions()->size(); i++) {
    ft_dims.push_back(i);
  }

  cuNDFFT<REAL>::instance()->ifft(&tmp, &ft_dims);

  if (!accumulate) 
    clear(out);
  
  this->mult_csm_conj_sum(&tmp,out);
}

//
// Instantiations
//

template class EXPORTGPUPMRI cuCartesianSenseOperator<float,1>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<float,2>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<float,3>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<float,4>;

template class EXPORTGPUPMRI cuCartesianSenseOperator<double,1>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<double,2>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<double,3>;
template class EXPORTGPUPMRI cuCartesianSenseOperator<double,4>;

