#pragma once

#include "cuCGMatrixOperator.h"
#include "ndarray_vector_td_utilities.h"

#include <cmath>
#include <boost/smart_ptr.hpp>

template <class REAL, class T> 
class cuCGImageOperator : public cuCGMatrixOperator<REAL,T>
{

 public:

  cuCGImageOperator( int device = -1 ) : cuCGMatrixOperator<REAL,T>( device ) {}
  virtual ~cuCGImageOperator() {}

  inline int set_encoding_operator( boost::shared_ptr< cuCGMatrixOperator<REAL,T> > encoding_operator )
  {
    encoding_operator_ = encoding_operator;
    return 0;
  }
  
  inline int compute( cuNDArray<T> *encoded_image, std::vector<unsigned int> *decoded_image_dimensions ) 
  { 
    if( !encoding_operator_.get() ){
      std::cout << std::endl << "cuCGImageOperator::compute: encoding operator not set" << std::endl;
      return -1;
    }
    
    // Reconstruct regularization image
    cuNDArray<T> tmp; tmp.create(decoded_image_dimensions, this->device_);
    encoding_operator_->mult_MH( encoded_image, &tmp );
    
    // Normalize to an average energy of "one intensity unit per image element"
    REAL sum = cuNDA_asum<REAL,T>( &tmp );
    REAL scale = ( (REAL) tmp.get_number_of_elements()/sum );
    cuNDA_scal<T>( scale*get_one<T>(), &tmp );

    // Square and reciprocalize image
    image_ = cuNDA_norm_squared<REAL,T>(&tmp, CUNDA_NDARRAY_DEVICE, CUNDA_NDARRAY_DEVICE);     
    cuNDA_reciprocal<REAL>(image_.get());
    
    return 0;
  }
    
  virtual int mult_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false)
  {
    std::cout << std::endl << "cuCGImageOperator::mult_M not defined." << std::endl;
    return -1;
  }
  
  virtual int mult_MH(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false)
  {
    std::cout << std::endl << "cuCGImageOperator::mult_MH not defined." << std::endl;
    return -1;
  }
  
  virtual int mult_MH_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false)
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
  
  cuNDArray<REAL>* get() { return image_.get(); }
  
private:
  boost::shared_ptr< cuCGMatrixOperator<REAL,T> > encoding_operator_;
  boost::shared_ptr< cuNDArray<REAL> > image_;
};
