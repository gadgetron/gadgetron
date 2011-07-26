#include "cuCGImageOperator.h"
#include "ndarray_vector_td_utilities.h"

#include <cmath>


template<class REAL, class T> 
int cuCGImageOperator<REAL,T>::compute( cuNDArray<T> *encoded_image, std::vector<unsigned int> *decoded_image_dimensions ) 
{ 
  if( !encoding_operator_.get() ){
    std::cout << std::endl << "cuCGImageOperator::compute: encoding operator not set" << std::endl;
    return -1;
  }
  
  // Reconstruct regularization image
  cuNDArray<T> tmp; tmp.create(decoded_image_dimensions, this->device_);
  if( encoding_operator_->mult_MH( encoded_image, &tmp ) != 0 ){
    std::cout << std::endl << "cuCGImageOperator::compute: adjoint encoding operator failed" << std::endl;
    return -1;
  }
  
  // Normalize to an average energy of "one intensity unit per image element"
  REAL sum = cuNDA_asum<REAL,T>( &tmp );
  REAL scale = ( (REAL) tmp.get_number_of_elements()/sum );
  cuNDA_scal<T>( scale*get_one<T>(), &tmp );
  
  // Square and reciprocalize image
  image_ = cuNDA_norm_squared<REAL,T>(&tmp, CUNDA_NDARRAY_DEVICE, CUNDA_NDARRAY_DEVICE);     
  cuNDA_reciprocal<REAL>(image_.get());
  
  return 0;
}

template<class REAL, class T> 
int cuCGImageOperator<REAL,T>::mult_MH_M( cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate )
{    
  int ret1 = this->set_device();
  bool ret2 = true;
  
  if( ret1 == 0 && !accumulate ) 
    ret2 = cuNDA_clear( out, get_zero<T>(), CUNDA_CURRENT_DEVICE );
  
  if( ret1 == 0 && ret2 )
    ret2 = cuNDA_axpy<REAL>( image_.get(), in, out, CUNDA_CURRENT_DEVICE );
  else 
    ret2 = false;
  
  ret1 = this->restore_device();
  
  if( ret1 == 0 && ret2 )
    return 0;
  else{
    std::cout << std::endl << "cuCGImageOperator::mult_MH_M failed" << std::endl;
    return -1;
  }
}

// Instantiation

template class cuCGImageOperator<float,float>;
template class cuCGImageOperator<float,float_complext::Type>;

template class cuCGImageOperator<double,double>;
template class cuCGImageOperator<double,double_complext::Type>;
