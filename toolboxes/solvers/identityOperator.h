#pragma once

#include "matrixOperator.h"
#include <iostream>

template <class REAL, class ARRAY_TYPE> class identityOperator 
	: public matrixOperator<REAL, ARRAY_TYPE>
{
 public:

  identityOperator() : matrixOperator<REAL,ARRAY_TYPE>() {}
  virtual ~identityOperator() {}
  
  virtual bool operator_xpy( ARRAY_TYPE*, ARRAY_TYPE* ) = 0;

  virtual int mult_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false)
  {
    // We will do only the most basic dimensionality checking
    if( in->get_number_of_elements() != out->get_number_of_elements() ){
      std::cout << std::endl << "Error :: identityOperator :: in/out dimensions mismatch" << std::endl;
      return -1;
    }

    bool ret = true;    
    
    if( accumulate )
      ret = operator_xpy( in, out );
    else 
      *out = *in;
    
    if( ret )      
      return 0;
    else{
      std::cout << std::endl << "Error :: identityOperator :: mult failed" << std::endl;
      return -1;
    }
  }
  
  virtual int mult_MH( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {
    return mult_M(in, out, accumulate);
  }
  
  virtual int mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out, bool accumulate = false )
  {
    return mult_M(in, out, accumulate);
  }
};
