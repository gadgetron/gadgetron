#include "cgOperatorSenseRHSBuffer.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template<class REAL, unsigned int D> int 
cgOperatorSenseRHSBuffer<REAL,D>::mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  std::cout << "cgOperatorSenseRHSBuffer::mult_M is not defined" << std::endl;
  return -1;
}

template<class REAL, unsigned int D> int 
cgOperatorSenseRHSBuffer<REAL,D>::mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate )
{
  if( !this->csm_.get() ){
    std::cerr << "cgOperatorSenseRHSBuffer::mult_MH: csm not set" << std::endl;
    return -1;
  }
  
  if( out->get_number_of_dimensions() != D && in->get_number_of_dimensions() != D+1 ){
    std::cerr << "cgOperatorSenseRHSBuffer::mult_MH array dimensions mismatch" << std::endl;
    return -1;
  }
  
  if( in->get_number_of_elements() != this->csm_->get_number_of_elements() ){
    std::cerr << "cgOperatorSenseRHSBuffer::mult_MH input array dimensions mismatch csm" << std::endl;
    return -2;
  }
  
  typename uintd<D>::Type image_size_in = vector_to_uintd<D>(in->get_dimensions());
  typename uintd<D>::Type image_size_out = vector_to_uintd<D>(out->get_dimensions());
  
  if( prod(image_size_in) != prod(image_size_out) ){
    std::cerr << "cgOperatorSenseRHSBuffer::mult_MH image dimensions mismatch" << std::endl;
    return -3;
  }

  if( !accumulate )
    cuNDA_clear<_complext>( out );
  
  if( mult_csm_conj_sum( in, out ) < 0 ) {
    std::cerr << "cgOperatorSenseRHSBuffer::mult_MH: Unable to multiply with conjugate of sensitivity maps and sum" << std::endl;
    return -4; 
  }
  
  return 0;
}

//
// Instantiations
//

template class cgOperatorSenseRHSBuffer<float,2>;
