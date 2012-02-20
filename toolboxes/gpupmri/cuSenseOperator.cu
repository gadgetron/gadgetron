#include "cuSenseOperator.h"
#include "vector_td_utilities.h"

template<class REAL> __global__ void 
mult_csm_kernel( complext<REAL> *in, complext<REAL> *out, complext<REAL> *csm,
		 unsigned int image_elements, unsigned int nframes, unsigned int ncoils )
{
  unsigned int idx = blockIdx.x*blockDim.x+threadIdx.x;
  if( idx < image_elements) {
    complext<REAL> _in = in[idx+blockIdx.y*image_elements];
    for( unsigned int i=0; i<ncoils; i++) {
      out[idx + blockIdx.y*image_elements + i*image_elements*nframes] =  _in * csm[idx+i*image_elements];
    }
  }
}

template<class REAL, unsigned int D> int 
cuSenseOperator<REAL,D>::mult_csm( cuNDArray<complext<REAL> >* in, cuNDArray<complext<REAL> >* out )
{  
  int device;
  if( cudaGetDevice( &device ) != cudaSuccess ){
    std::cerr << "cuOperatorSense::mult_csm: unable to query current device" << std::endl;
    return -1;
  }
  
  if( !in || in->get_device() != device || !out || out->get_device() != device || !this->csm_.get() || this->csm_->get_device() != device ){
    std::cerr << "cuOperatorSense::mult_csm: input array not residing current device" << std::endl;
    return -1;
  }
  
  if( in->get_number_of_dimensions() < D  || in->get_number_of_dimensions() > D+1 ){
    std::cerr << "cuOperatorSense::mult_csm: unexpected input dimensionality" << std::endl;
    return -1;
  }

  if( in->get_number_of_dimensions() > out->get_number_of_dimensions() ){
    std::cerr << "cuOperatorSense::mult_csm: input dimensionality cannot exceed output dimensionality" << std::endl;
    return -1;
  }

  if( this->csm_->get_number_of_dimensions() != D+1 ) {
    std::cerr << "cuOperatorSense::mult_csm: input dimensionality of csm not as expected" << std::endl;
    return -1;
  }

  unsigned int num_image_elements = 1;
  for( unsigned int d=0; d<D; d++ )
    num_image_elements *= in->get_size(d);
  
  unsigned int num_frames = in->get_number_of_elements() / num_image_elements;
  
  dim3 blockDim(256);
  dim3 gridDim((num_image_elements+blockDim.x-1)/blockDim.x, num_frames);

  mult_csm_kernel<REAL><<< gridDim, blockDim >>>
    ( in->get_data_ptr(), out->get_data_ptr(), this->csm_->get_data_ptr(), num_image_elements, num_frames, this->csm_->get_size(D) );

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cuOperatorSense::mult_csm: unable to multiply with coil sensitivities: " << 
      cudaGetErrorString(err) << std::endl;
    return -2;
  }

  return 0;
}

template <class REAL> __global__ void 
mult_csm_conj_sum_kernel( complext<REAL> *in, complext<REAL> *out, complext<REAL> *csm,
						  unsigned int image_elements, unsigned int nframes, unsigned int ncoils )
{
  unsigned int idx = blockIdx.x*blockDim.x+threadIdx.x;
  if( idx < image_elements ) {
    complext<REAL> _out =complext<REAL>(0);
    for( unsigned int i = 0; i < ncoils; i++ ) {
      _out += in[idx+blockIdx.y*image_elements+i*nframes*image_elements] * conj(csm[idx+i*image_elements]);
    }
    out[idx+blockIdx.y*image_elements] = _out;
  }
}

template<class REAL, unsigned int D> int 
cuSenseOperator<REAL,D>::mult_csm_conj_sum( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out )
{
  int device;
  if( cudaGetDevice( &device ) != cudaSuccess ){
    std::cerr << "cuOperatorSense::mult_csm_conj_sum: unable to query current device" << std::endl;
    return -1;
  }
  
  if( !in || in->get_device() != device || !out || out->get_device() != device || !this->csm_.get() || this->csm_->get_device() != device ){
    std::cerr << "cuOperatorSense::mult_csm_conj_sum: input array not residing current device" << std::endl;
    return -1;
  }
  
  if( out->get_number_of_dimensions() < D  || out->get_number_of_dimensions() > D+1 ){
    std::cerr << "cuOperatorSense::mult_csm_conj_sum: unexpected output dimensionality" << std::endl;
    return -1;
  }

  if( out->get_number_of_dimensions() > in->get_number_of_dimensions() ){
    std::cerr << "cuOperatorSense::mult_csm_conj_sum: output dimensionality cannot exceed input dimensionality" << std::endl;
    return -1;
  }

  if( this->csm_->get_number_of_dimensions() != D+1 ) {
    std::cerr << "cuOperatorSense::mult_csm_conj_sum: input dimensionality of csm not as expected" << std::endl;
    return -1;
  }

  unsigned int num_image_elements = 1;
  for( unsigned int d=0; d<D; d++ )
    num_image_elements *= out->get_size(d);
  
  unsigned int num_frames = out->get_number_of_elements() / num_image_elements;

  dim3 blockDim(256);
  dim3 gridDim((num_image_elements+blockDim.x-1)/blockDim.x, num_frames);

  mult_csm_conj_sum_kernel<REAL><<< gridDim, blockDim >>>
    ( in->get_data_ptr(), out->get_data_ptr(), this->csm_->get_data_ptr(), num_image_elements, num_frames, this->csm_->get_size(D) );

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cuOperatorSense::mult_csm_conj_sum: unable to combine coils " << 
      cudaGetErrorString(err) << std::endl;
    return -2;
  }
  
  return 0;
}


//
// Instantiations
//

template class EXPORTGPUPMRI cuSenseOperator<float,2>;
template class EXPORTGPUPMRI cuSenseOperator<float,3>;
template class EXPORTGPUPMRI cuSenseOperator<float,4>;

template class EXPORTGPUPMRI cuSenseOperator<double,2>;
template class EXPORTGPUPMRI cuSenseOperator<double,3>;
template class EXPORTGPUPMRI cuSenseOperator<double,4>;
