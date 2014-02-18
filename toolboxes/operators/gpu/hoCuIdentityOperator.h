/** \file hoCuIdentityOperator.h
    \brief Instantiation of the identity operator for array type hoCuNDArray
    
    The file hoCuIdentityOperator.h is a convienience wrapper for the device independent identityOperator class.
    The class hoCuIdentityOperator instantiates the identityOperator for the hoCuNDArray
    and the header furthermore includes additional neccessary header files.
*/

#pragma once

#include "hoCuNDArray_math.h"
#include "identityOperator.h"

namespace Gadgetron{
  
  /** \class hoCuIdentityOperator
      \brief Instantiation of the identity operator for array type hoCuNDArray
      
      The class hoCuIdentityOperator is a convienience wrapper for the device independent identityOperator.
      hoCuIdentityOperator instantiates the identityOperator for type hoCuNDArray<T>.
  */
  template <class T> class hoCuIdentityOperator : public identityOperator< hoCuNDArray<T> >
  {
  public:    
    hoCuIdentityOperator() : identityOperator< hoCuNDArray<T> >() {}
    virtual ~hoCuIdentityOperator() {}
  }; 
}
