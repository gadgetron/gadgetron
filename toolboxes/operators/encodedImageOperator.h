/** \file encodedImageOperator.h
    \brief Regularization operator for encoded images. Careful, only implements mult_MH_M and not (yet) mult_M and mult_MH.
*/

#pragma once

#include "imageOperator.h"

namespace Gadgetron{

  template <class ARRAY_TYPE_REAL, class ARRAY_TYPE_OPERATOR> class encodedImageOperator
    : public imageOperator<ARRAY_TYPE_REAL, ARRAY_TYPE_OPERATOR>
  {
  
  public:
  
    encodedImageOperator() : imageOperator<ARRAY_TYPE_REAL, ARRAY_TYPE_OPERATOR>() {}
    virtual ~encodedImageOperator() {}
 
    // Set encoding operator for the regularization image
    virtual void set_encoding_operator( boost::shared_ptr< linearOperator<ARRAY_TYPE_OPERATOR> > encoding_operator )

    {
      encoding_operator_ = encoding_operator;
    }
  
    // Apply regularization image operator
    virtual void mult_MH_M( ARRAY_TYPE_OPERATOR *in, ARRAY_TYPE_OPERATOR *out, bool accumulate = false )
    {    
      if( !encoding_operator_.get() ){
        throw std::runtime_error("encodedImageOperator::mult_MH_M failed : encoding operator not set");
      }
    
      ARRAY_TYPE_OPERATOR tmp(in->get_dimensions());

      encoding_operator_->mult_M( in, &tmp );
 
      ARRAY_TYPE_OPERATOR tmp2(in->get_dimensions());

      imageOperator<ARRAY_TYPE_REAL, ARRAY_TYPE_OPERATOR>::mult_MH_M( &tmp, &tmp2 );
    
      encoding_operator_->mult_MH( &tmp2, out, accumulate );
    }  
  
  private:
    boost::shared_ptr< linearOperator<ARRAY_TYPE_OPERATOR> > encoding_operator_;
  };
}
