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

  if( out->get_number_of_dimensions() < D || out->get_number_of_dimensions() > D+1 ){
    std::cerr << "cgOperatorSenseRHSBuffer::mult_MH: unexpected output array dimensions" << std::endl;
    return -1;
  }

  if( out->get_number_of_dimensions() > in->get_number_of_dimensions() ){
    std::cerr << "cgOperatorSenseRHSBuffer::mult_MH: array dimensions mismatch" << std::endl;
    return -1;
  }
    
  for( unsigned int d=0; d<out->get_number_of_dimensions(); d++ ){
    if( out->get_size(d) != in->get_size(d) ){
      std::cerr << "cgOperatorSenseRHSBuffer::mult_MH: array dimensions mismatch" << std::endl;
      return -3;
    }
    if( d<D && in->get_size(d) != this->csm_->get_size(d) ){
      std::cerr << "cgOperatorSenseRHSBuffer::mult_MH: input array dimensions mismatch the csm" << std::endl;
      return -3;
    }
  }
  
  if( in->get_size(in->get_number_of_dimensions()-1) != 
      this->csm_->get_size(this->csm_->get_number_of_dimensions()-1 )){
    std::cerr << "cgOperatorSenseRHSBuffer::mult_MH: the number of coils cannot differ between the input array and the csm" << std::endl;
    return -1;
  }
  
  if( !accumulate ){
    int ret1 = this->set_device();
    if( ret1 == 0 ) cuNDA_clear<_complext>( out );
    int ret2 = this->restore_device();
    
    if( ret1<0 || ret2<0 ){
      std::cerr << "cgOperatorSenseRHSBuffer::mult_MH: device error" << std::endl;
      return -4; 
    }
  }
 
  if( mult_csm_conj_sum( in, out ) < 0 ) {
    std::cerr << "cgOperatorSenseRHSBuffer::mult_MH: Unable to multiply with conjugate of sensitivity maps and sum" << std::endl;
    return -4; 
  }
  
  return 0;
}

//
// Instantiations
//

template class EXPORTGPUPMRI cgOperatorSenseRHSBuffer<float,2>;
