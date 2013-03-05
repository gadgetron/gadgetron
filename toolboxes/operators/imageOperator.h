#pragma once

#include "linearOperator.h"
#include "vector_td_utilities.h"

#include <boost/smart_ptr.hpp>
#include <vector>

namespace Gadgetron{


template <class ARRAY_TYPE_REAL, class ARRAY_TYPE_OPERATOR> class imageOperator
  : public linearOperator<ARRAY_TYPE_OPERATOR>
{
  typedef typename ARRAY_TYPE_REAL::element_type REAL;
  typedef typename ARRAY_TYPE_OPERATOR::element_type ELEMENT_TYPE;
public:
  
  imageOperator() : linearOperator<ARRAY_TYPE_OPERATOR>() {}
  virtual ~imageOperator() {}
  
  // Get regularization image
  virtual ARRAY_TYPE_REAL* get() { return image_.get(); }
    
  // Compute regularization image (apply the adjoint encoding operator on the encoded image)
  virtual void compute( ARRAY_TYPE_OPERATOR *image )
  {     
    // Make temporary copy of input
    ARRAY_TYPE_OPERATOR tmp(*image);

    // Normalize to an average energy of "one intensity unit per image element"
    REAL sum = asum( &tmp );
    REAL scale = ( (REAL) tmp.get_number_of_elements()/sum );
    tmp *= scale;

    
    // Reciprocalize image
    image_ =  abs(&tmp);
    image_->reciprocal();

  }
  
  // Apply regularization image operator
  virtual void mult_MH_M( ARRAY_TYPE_OPERATOR *in, ARRAY_TYPE_OPERATOR *out, bool accumulate = false )
  {    
    
  	ARRAY_TYPE_OPERATOR * tmp;
    if( !accumulate ){
    	tmp = out;
    	*tmp = *in;
    } else
    	tmp = new ARRAY_TYPE_OPERATOR(*in);
    
    *tmp *= *image_;
    *tmp *= ELEMENT_TYPE(2);

    if (accumulate){
    	*out += *tmp;
    	delete tmp;
    }
  }
  
  virtual void mult_M( ARRAY_TYPE_OPERATOR *in, ARRAY_TYPE_OPERATOR *out, bool accumulate = false )
  {

  	ARRAY_TYPE_OPERATOR * tmp;
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
  
  virtual boost::shared_ptr< linearOperator<ARRAY_TYPE_OPERATOR > > clone()
   {
     return linearOperator<ARRAY_TYPE_OPERATOR>::clone(this);
   }

protected:
  boost::shared_ptr< ARRAY_TYPE_REAL > image_;
};
}
