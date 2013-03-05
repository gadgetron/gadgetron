#pragma once

#include "linearOperator.h"
#include <iostream>

namespace Gadgetron{
template < class ARRAY_TYPE> class identityOperator
	: public linearOperator<ARRAY_TYPE>
{
 public:

  identityOperator() : linearOperator<ARRAY_TYPE>() {}
  virtual ~identityOperator() {}
  


  virtual void mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {
    // We will do only the most basic dimensionality checking
    if( in->get_number_of_elements() != out->get_number_of_elements() ){
      BOOST_THROW_EXCEPTION(runtime_error("Error :: identityOperator :: in/out dimensions mismatch"));
    }


    
    if( accumulate )
    	*out += *in;

    else 
      *out = *in;
    

  }
  
  virtual void mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {
    return mult_M(in, out, accumulate);
  }
  
  virtual void mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {
    return mult_M(in, out, accumulate);
  }

  virtual boost::shared_ptr< linearOperator< ARRAY_TYPE > > clone()
   {
     return linearOperator<  ARRAY_TYPE >::clone(this);
   }
};
}
