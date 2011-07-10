#pragma once

#include "cuCGMatrixOperator.h"
#include "gadgetron_export.h"

#include <boost/smart_ptr.hpp>

template <class REAL, class T> 
class EXPORTGPUCG cuCGImageOperator : public cuCGMatrixOperator<REAL,T>
{

 public:

  cuCGImageOperator( int device = -1 ) : cuCGMatrixOperator<REAL,T>( device ) {}
  virtual ~cuCGImageOperator() {}
  
  // Set encoding operator for the regularization image
  inline void set_encoding_operator( boost::shared_ptr< cuCGMatrixOperator<REAL,T> > encoding_operator ) {
    encoding_operator_ = encoding_operator;
  }
  
  // Compute regularization image (apply the adjoint encoding operator on the encoded image)
  virtual int compute( cuNDArray<T> *encoded_image, std::vector<unsigned int> *decoded_image_dimensions );
  
  // Get regularization image
  inline cuNDArray<REAL>* get() { return image_.get(); }
  
  // Apply regularization image operator
  virtual int mult_MH_M( cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false );
 
  virtual int mult_M( cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false ){
    std::cout << std::endl << "cuCGImageOperator::mult_M not defined." << std::endl;
    return -1;
  }
  
  virtual int mult_MH( cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false ){
    std::cout << std::endl << "cuCGImageOperator::mult_MH not defined." << std::endl;
    return -1;
  }
  
private:
  boost::shared_ptr< cuCGMatrixOperator<REAL,T> > encoding_operator_;
  boost::shared_ptr< cuNDArray<REAL> > image_;
};
