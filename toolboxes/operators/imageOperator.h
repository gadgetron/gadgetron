/** \file imageOperator.h
    \brief Base class for the image regularization operators.
*/

#pragma once

#include "linearOperator.h"
#include "GadgetronTimer.h"

namespace Gadgetron{
  
  template <class ARRAY_TYPE_REAL, class ARRAY_TYPE_OPERATOR> class imageOperator : public linearOperator<ARRAY_TYPE_OPERATOR>
  {

  protected:
    typedef typename ARRAY_TYPE_REAL::element_type REAL;
    typedef typename ARRAY_TYPE_OPERATOR::element_type ELEMENT_TYPE;
    
  public:
    
    imageOperator() : linearOperator<ARRAY_TYPE_OPERATOR>(), offset_(REAL(0)) {}
    virtual ~imageOperator() {}
  
    // Get regularization image
    virtual boost::shared_ptr<ARRAY_TYPE_REAL> get() { return image_; }
    
    // Compute regularization image
    virtual void compute( ARRAY_TYPE_OPERATOR *image, bool offset_estimation = true )
    {
      // Make temporary copy of input
      ARRAY_TYPE_OPERATOR tmp(*image);

      // Normalize to an average energy of "one intensity unit per image element"
      REAL sum = asum( &tmp );
      REAL scale = ( (REAL) tmp.get_number_of_elements()/sum );
      tmp *= scale;

      image_ =  abs(&tmp);

      if( offset_estimation )
	offset_ = estimate_offset();
      
      // Reciprocalize image
      if(offset_ > REAL(0)) *image_ += offset_;      
      reciprocal_inplace(image_.get());
    }
    
    // Apply regularization image operator
    virtual void mult_MH_M( ARRAY_TYPE_OPERATOR *in, ARRAY_TYPE_OPERATOR *out, bool accumulate = false )
    {        
      ARRAY_TYPE_OPERATOR *tmp;
      if( !accumulate ){
    	tmp = out;
    	*tmp = *in;
      } 
      else
    	tmp = new ARRAY_TYPE_OPERATOR(*in);
      
      *tmp *= *image_;
      *tmp *= *image_;
      
      if (accumulate){
    	*out += *tmp;
    	delete tmp;
      }
    }
  
    virtual void mult_M( ARRAY_TYPE_OPERATOR *in, ARRAY_TYPE_OPERATOR *out, bool accumulate = false )
    {
      ARRAY_TYPE_OPERATOR *tmp;
      if( !accumulate ){
	tmp = out;
	*tmp = *in;
      } else
	tmp = new ARRAY_TYPE_OPERATOR(*in);

      *tmp *= *image_;

      if (accumulate){
	*out += *tmp;
	delete tmp;
      }
    }
  
    virtual void mult_MH( ARRAY_TYPE_OPERATOR *in, ARRAY_TYPE_OPERATOR *out, bool accumulate = false )
    {
      mult_M(in,out,accumulate);
    }

  
  protected:
    // Estimate offset to the regularization image
    virtual REAL estimate_offset()=0;

  protected:
    boost::shared_ptr< ARRAY_TYPE_REAL > image_;
    REAL offset_;
  };
}
