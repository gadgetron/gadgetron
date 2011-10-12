#pragma once

#include "imageOperator.h"

template <class REAL, class ARRAY_TYPE_REAL, class ARRAY_TYPE_OPERATOR> class encodedImageOperator : public imageOperator<REAL, ARRAY_TYPE_REAL, ARRAY_TYPE_OPERATOR>
{
  
public:
  
  encodedImageOperator() : imageOperator<REAL, ARRAY_TYPE_REAL, ARRAY_TYPE_OPERATOR>() {}
  virtual ~encodedImageOperator() {}
 
  // Set encoding operator for the regularization image
  inline void set_encoding_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_OPERATOR> > encoding_operator )				     

  {
    encoding_operator_ = encoding_operator;
  }
  
  // Apply regularization image operator
  virtual int mult_MH_M( ARRAY_TYPE_OPERATOR *in, ARRAY_TYPE_OPERATOR *out, bool accumulate = false )
  {    
    if( !encoding_operator_.get() ){
      std::cout << std::endl << "encodedImageOperator::mult_MH_M failed : encoding operator not set" << std::endl;
      return -1;	
    }
    
    ARRAY_TYPE_OPERATOR tmp; 
    if( tmp.create( in->get_dimensions().get() ) < 0 ) {
      std::cout << std::endl << "encodedImageOperator::mult_MH_M : decoded image allocation failed." << std::endl;
      return -1;	
    }

    if( encoding_operator_->mult_M( in, &tmp ) < 0 ){
      std::cout << std::endl << "encodedImageOperator::mult_MH_M : forwards encoding operator failed." << std::endl;
      return -1;	
    }
    ARRAY_TYPE_OPERATOR tmp2; 
    if( tmp2.create( in->get_dimensions().get() ) < 0 ) {
      std::cout << std::endl << "encodedImageOperator::mult_MH_M : decoded image allocation failed." << std::endl;
      return -1;	
    }
    if( imageOperator<REAL, ARRAY_TYPE_REAL, ARRAY_TYPE_OPERATOR>::mult_MH_M( &tmp, &tmp2 ) < 0 ){
      std::cout << std::endl << "encodedImageOperator::mult_MH_M : error from inherited class" << std::endl;
      return -1;	
    }
    
    if( encoding_operator_->mult_MH( &tmp2, out, accumulate ) < 0 ){
      std::cout << std::endl << "encodedImageOperator::mult_MH_M : backwards encoding operator failed." << std::endl;
      return -1;	
    }

    return 0;
  }  
  
private:
  boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_OPERATOR> > encoding_operator_;
};
