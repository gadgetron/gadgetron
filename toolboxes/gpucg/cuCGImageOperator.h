#pragma once

#include "cuCGMatrixOperator.h"
#include "ndarray_vector_td_utilities.h"

#include <cublas_v2.h>
#include <cmath>

template <class REAL, class T> 
class cuCGImageOperator : public cuCGMatrixOperator<T>
{
 public:

  cuCGImageOperator( cuNDArray<T> *encoded_image, const std::vector<unsigned int> &decoded_image_dimensions,
		     cuCGMatrixOperator<T> *encoding_matrix, cublasHandle_t handle ) 
  { 
    // Reconstruct regularization image
    cuNDArray<T> tmp; tmp.create(decoded_image_dimensions);
    encoding_matrix->mult_MH( encoded_image, &tmp );
    
    // Normalize to an average energy of "one intensity unit per image element"
    REAL sum = cuNDA_asum<REAL,T>( &tmp, handle );
    REAL scale = ( (REAL) tmp.get_number_of_elements()/sum );
    cuNDA_scal<T>( scale*get_one<T>(), &tmp, handle );

    // Square and reciprocalize image
    image_ = cuNDA_norm_squared<REAL,T>(&tmp);     
    cuNDA_reciprocal<REAL>(image_.get()); 
  }

  virtual ~cuCGImageOperator() {}
  
  virtual int mult_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false){
    std::cout << std::endl << "cuCGImageOperator::mult_M not defined." << std::endl;
    return -1;
  }
  
  virtual int mult_MH(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false){
    std::cout << std::endl << "cuCGImageOperator::mult_MH not defined." << std::endl;
    return -1;
  }
  
  virtual int mult_MH_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false){

    if( !accumulate ) 
      cuNDA_clear(out);
    
    if( cuNDA_axpy<REAL>( image_.get(), in, out ) )
      return 0;
    else
      return -2;
  }
  
  cuNDArray<REAL>* get() { return image_.get(); }
  
private:
  std::auto_ptr< cuNDArray<REAL> > image_;
};
