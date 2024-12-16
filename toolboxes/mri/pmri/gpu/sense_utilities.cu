#include "sense_utilities.h"
#include "vector_td_utilities.h"
#include <sstream>

namespace Gadgetron{

  template<class REAL> __global__ void
  mult_csm_kernel( const complext<REAL> * __restrict__ in, complext<REAL> * __restrict__ out, complext<REAL> *csm,
		   size_t image_elements, unsigned int nframes, unsigned int ncoils )
  {
    unsigned int idx = blockIdx.x*blockDim.x+threadIdx.x;
    if( idx < image_elements) {
      complext<REAL> _in = in[idx+blockIdx.y*image_elements];
      for( unsigned int i=0; i<ncoils; i++) {
	out[idx + blockIdx.y*image_elements + i*image_elements*nframes] =  _in * csm[idx+i*image_elements];
      }
    }
  }

  template<class REAL, unsigned int D> void
  csm_mult_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, cuNDArray< complext<REAL> > *csm )
  {
    int device;
    if( cudaGetDevice( &device ) != cudaSuccess ){
      throw cuda_error( "mult_csm: unable to query current device");
    }

    if( !in || in->get_device() != device || !out || out->get_device() != device || !csm || csm->get_device() != device ){
      throw cuda_error("mult_csm: array not residing current device");
    }

    if( in->get_number_of_dimensions() < D  || in->get_number_of_dimensions() > D+1 ){
      throw std::runtime_error("mult_csm: unexpected input dimensionality");
    }

    if( in->get_number_of_dimensions() > out->get_number_of_dimensions() ){
      throw std::runtime_error("mult_csm: input dimensionality cannot exceed output dimensionality");
    }

    if( csm->get_number_of_dimensions() != D+1 ) {
      throw std::runtime_error("mult_csm: input dimensionality of csm not as expected");
    }

    unsigned int num_image_elements = 1;
    for( unsigned int d=0; d<D; d++ )
      num_image_elements *= in->get_size(d);

    unsigned int num_frames = in->get_number_of_elements() / num_image_elements;

    dim3 blockDim(256);
    dim3 gridDim((num_image_elements+blockDim.x-1)/blockDim.x, num_frames);

    mult_csm_kernel<REAL><<< gridDim, blockDim >>>
      ( in->get_data_ptr(), out->get_data_ptr(), csm->get_data_ptr(), num_image_elements, num_frames, csm->get_size(D) );

    cudaError_t err = cudaGetLastError();
    if( err != cudaSuccess ){
      std::stringstream ss;
      ss << "mult_csm: unable to multiply with coil sensitivities: " <<
	cudaGetErrorString(err);
      throw cuda_error(ss.str());

    }
  }

  template <class REAL> __global__ void
  mult_csm_conj_sum_kernel(const  complext<REAL> * __restrict__ in, complext<REAL> * __restrict__ out, const complext<REAL> * __restrict__ csm,
			    size_t image_elements, unsigned int nframes, unsigned int ncoils )
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

  template<class REAL, unsigned int D> void
  csm_mult_MH( cuNDArray<complext<REAL> > *in, cuNDArray<complext<REAL> > *out, cuNDArray<complext<REAL> > *csm )
  {
    int device;
    if( cudaGetDevice( &device ) != cudaSuccess ){
      throw cuda_error("mult_csm_conj_sum: unable to query current device");
    }

    if( !in || in->get_device() != device || !out || out->get_device() != device || !csm || csm->get_device() != device ){
      throw std::runtime_error("mult_csm_conj_sum: array not residing current device");
    }

    if( out->get_number_of_dimensions() < D  || out->get_number_of_dimensions() > D+1 ){
      throw std::runtime_error("mult_csm_conj_sum: unexpected output dimensionality");
    }

    if( out->get_number_of_dimensions() > in->get_number_of_dimensions() ){
      throw std::runtime_error("mult_csm_conj_sum: output dimensionality cannot exceed input dimensionality");
    }

    if( csm->get_number_of_dimensions() != D+1 ) {
      throw std::runtime_error("mult_csm_conj_sum: input dimensionality of csm not as expected");
    }

    unsigned int num_image_elements = 1;
    for( unsigned int d=0; d<D; d++ )
      num_image_elements *= out->get_size(d);

    unsigned int num_frames = out->get_number_of_elements() / num_image_elements;

    dim3 blockDim(256);
    dim3 gridDim((num_image_elements+blockDim.x-1)/blockDim.x, num_frames);

    mult_csm_conj_sum_kernel<REAL><<< gridDim, blockDim >>>
      ( in->get_data_ptr(), out->get_data_ptr(), csm->get_data_ptr(), num_image_elements, num_frames, csm->get_size(D) );

    cudaError_t err = cudaGetLastError();
    if( err != cudaSuccess ){
      std::stringstream ss;
      ss << "mult_csm_conj_sum: unable to combine coils " <<
	cudaGetErrorString(err);
      throw cuda_error(ss.str());
    }
  }

  // Instantiation

  template void csm_mult_M<float,1>( cuNDArray< complext<float> >*, cuNDArray< complext<float> >*, cuNDArray< complext<float> >*);
  template void csm_mult_M<float,2>( cuNDArray< complext<float> >*, cuNDArray< complext<float> >*, cuNDArray< complext<float> >*);
  template void csm_mult_M<float,3>( cuNDArray< complext<float> >*, cuNDArray< complext<float> >*, cuNDArray< complext<float> >*);
  template void csm_mult_M<float,4>( cuNDArray< complext<float> >*, cuNDArray< complext<float> >*, cuNDArray< complext<float> >*);

  template void csm_mult_M<double,1>( cuNDArray< complext<double> >*, cuNDArray< complext<double> >*, cuNDArray< complext<double> >*);
  template void csm_mult_M<double,2>( cuNDArray< complext<double> >*, cuNDArray< complext<double> >*, cuNDArray< complext<double> >*);
  template void csm_mult_M<double,3>( cuNDArray< complext<double> >*, cuNDArray< complext<double> >*, cuNDArray< complext<double> >*);
  template void csm_mult_M<double,4>( cuNDArray< complext<double> >*, cuNDArray< complext<double> >*, cuNDArray< complext<double> >*);

  template void csm_mult_MH<float,1>( cuNDArray< complext<float> >*, cuNDArray< complext<float> >*, cuNDArray< complext<float> >*);
  template void csm_mult_MH<float,2>( cuNDArray< complext<float> >*, cuNDArray< complext<float> >*, cuNDArray< complext<float> >*);
  template void csm_mult_MH<float,3>( cuNDArray< complext<float> >*, cuNDArray< complext<float> >*, cuNDArray< complext<float> >*);
  template void csm_mult_MH<float,4>( cuNDArray< complext<float> >*, cuNDArray< complext<float> >*, cuNDArray< complext<float> >*);

  template void csm_mult_MH<double,1>( cuNDArray< complext<double> >*, cuNDArray< complext<double> >*, cuNDArray< complext<double> >*);
  template void csm_mult_MH<double,2>( cuNDArray< complext<double> >*, cuNDArray< complext<double> >*, cuNDArray< complext<double> >*);
  template void csm_mult_MH<double,3>( cuNDArray< complext<double> >*, cuNDArray< complext<double> >*, cuNDArray< complext<double> >*);
  template void csm_mult_MH<double,4>( cuNDArray< complext<double> >*, cuNDArray< complext<double> >*, cuNDArray< complext<double> >*);
}
