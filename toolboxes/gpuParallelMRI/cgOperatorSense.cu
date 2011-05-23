#include "cgOperatorSense.h"
#include "vector_td_utilities.h"

template<class REAL, unsigned int D> int 
cgOperatorSense<REAL,D>::set_csm( cuNDArray<_complext>* csm ) 
{
  if( csm && csm->get_number_of_dimensions() == D+1 ) {
    csm_ = csm;
    ncoils_ = csm_->get_size(D);
    dimensionsI_ = csm->get_dimensions();
    dimensionsI_.pop_back();
    return 0;
  }
  else{
    std::cerr << "cgOperatorSense::set_csm: dimensionality mismatch " << std::endl;
    return -1;
  }
}

template<class REAL> __global__ void 
mult_csm_kernel( typename complext<REAL>::Type* in, typename complext<REAL>::Type* out, typename complext<REAL>::Type* csm, 
		 unsigned long image_elements, unsigned int nframes, unsigned int ncoils )
{
  unsigned long idx = blockIdx.x*blockDim.x+threadIdx.x;
  if( idx < image_elements) {
    typedef typename complext<REAL>::Type T;
    T _in = in[idx+blockIdx.y*image_elements];
    for( unsigned int i=0; i<ncoils; i++) {
      out[idx + blockIdx.y*image_elements + i*image_elements*nframes] = mul<T>( _in, csm[idx+i*image_elements] );
    }
  }
}

template<class REAL, unsigned int D> int 
cgOperatorSense<REAL,D>::mult_csm( cuNDArray<_complext>* in, cuNDArray<_complext>* out )
{  
  // protected method: we skip dimensionality checks and trust the caller

  unsigned int num_image_elements = 1;
  for( unsigned int d=0; d<D; d++ )
    num_image_elements *= in->get_size(d);

  unsigned int num_frames = in->get_number_of_elements() / num_image_elements;

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)num_image_elements/blockDim.x), num_frames);
  mult_csm_kernel<REAL><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), csm_->get_data_ptr(), 
						  num_image_elements, num_frames, ncoils_ );
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cgOperatorSense::mult_csm Unable to multiply with coil sensitivities: " << 
      cudaGetErrorString(err) << std::endl;
    return -1;
  }
  return 0;
}

template <class REAL> __global__ void 
mult_csm_conj_sum_kernel( typename complext<REAL>::Type *in, typename complext<REAL>::Type *out, typename complext<REAL>::Type *csm, 
			  unsigned int image_elements, unsigned int nframes, unsigned int ncoils )
{
  unsigned int idx = blockIdx.x*blockDim.x+threadIdx.x;
  if( idx < image_elements ) {
      typedef typename complext<REAL>::Type T;
      T _out = get_zero<T>();
      for( unsigned int i = 0; i < ncoils; i++ ) {
	_out += mul<T>( in[idx+blockIdx.y*image_elements+i*nframes*image_elements], conj<T>(csm[idx+i*image_elements]) );
      }
      out[idx+blockIdx.y*image_elements] = _out;
  }
}

template<class REAL, unsigned int D> int 
cgOperatorSense<REAL,D>:: mult_csm_conj_sum( cuNDArray<_complext>* in, cuNDArray<_complext>* out )
{
  // protected method: we skip dimensionality checks and trust the caller

  unsigned int num_image_elements = 1;
  for( unsigned int d=0; d<D; d++ )
    num_image_elements *= out->get_size(d);

  unsigned int num_frames = out->get_number_of_elements() / num_image_elements;

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)num_image_elements/blockDim.x), num_frames);

  mult_csm_conj_sum_kernel<REAL><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), csm_->get_data_ptr(),
							   num_image_elements, num_frames, ncoils_ );
  
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cgOperatorSense::mult_EM : Unable to combine coils " << 
      cudaGetErrorString(err) << std::endl;
    return -1;
  }

  return 0;
}

template<class REAL, unsigned int D> int 
cgOperatorSense<REAL,D>::mult_MH_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  // Leave it to the inherited classes to validate the input

  cuNDArray<_complext> tmp;
  if( !tmp.create(dimensionsK_) ) {
    std::cerr << "cgOperatorSense::mult_MH_M: Unable to create temporary storage" << std::endl;
    return -1;
  }

  if( mult_M( in, &tmp, false ) < 0 ) {
    std::cerr << "cgOperatorSense::mult_MH_M: Unable to perform mult_M" << std::endl;
    return -2;
  }

  if( mult_MH( &tmp, out, accumulate ) < 0 ) {
    std::cerr << "cgOperatorSense::mult_MH_M: Unable to perform mult_MH" << std::endl;
    return -3;
  }

  return 0;
}

//
// Instantiations
//

template class cgOperatorSense<float,2>;
