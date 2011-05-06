#include "cgOperatorSense.h"
#include "vector_td_utilities.h"

template<class REAL, unsigned int D> int 
cgOperatorSense<REAL,D>::set_csm( cuNDArray<_complext>* csm ) 
{
  if( csm && csm->get_number_of_dimensions() == D+1 ) {
    csm_ = csm;
    ncoils_ = csm_->get_size(D);
    dimensions_ = csm->get_dimensions();
    dimensions_.pop_back();
    return 0;
  }
  else{
    std::cerr << "cgOperatorSense::set_csm: dimensionality mismatch " << std::endl;
    return -1;
  }
}

template<class REAL> __global__ void 
mult_csm_kernel( typename complext<REAL>::Type* in, typename complext<REAL>::Type* out, typename complext<REAL>::Type* csm, 
		 unsigned long image_elements, unsigned int ncoils )
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if( idx_in < image_elements) {
    for( unsigned int i=0; i<ncoils; i++) {
      typedef typename complext<REAL>::Type T;
      out[idx_in + i*image_elements] = mul<T,T>( in[idx_in], csm[idx_in+i*image_elements] );
    }
  }
}

template<class REAL, unsigned int D> int 
cgOperatorSense<REAL,D>::mult_csm( cuNDArray<_complext>* in, cuNDArray<_complext>* out )
{  
  // protected method: we skip dimensionality chekcs and trust the caller

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)in->get_number_of_elements()/blockDim.x));
  mult_csm_kernel<REAL><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), csm_->get_data_ptr(), 
						  in->get_number_of_elements(), ncoils_ );
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
			  unsigned int image_elements, unsigned int ncoils )
{
  unsigned int idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if( idx_in < image_elements ) {
    for( unsigned int i = 0; i < ncoils; i++ ) {
      typedef typename complext<REAL>::Type T;
      T tmp = mul<T,T>( in[idx_in+i*image_elements], conj<REAL>(csm[idx_in+i*image_elements]) );
      out[idx_in] += tmp;
    }
  }
}

template<class REAL, unsigned int D> int 
cgOperatorSense<REAL,D>:: mult_csm_conj_sum( cuNDArray<_complext>* in, cuNDArray<_complext>* out )
{
  // protected method: we skip dimensionality chekcs and trust the caller

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)out->get_number_of_elements()/blockDim.x));
  mult_csm_conj_sum_kernel<REAL><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), csm_->get_data_ptr(),
							   out->get_number_of_elements(), ncoils_ );
  
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
  if( !tmp.create(dimensions_out_) ) {
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
