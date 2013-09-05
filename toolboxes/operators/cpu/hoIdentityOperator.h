/** \file hoIdentityOperator.h
    \brief Instantiation of the identity operator on the cpu.
    
    The file hoIdentityOperator.h is a convienience wrapper for the device independent identityOperator class.
    The class hoIdentityOperator instantiates the identityOperator for the hoNDArray
    and the header furthermore includes additional neccessary header files.
*/

#pragma once

#include "hoNDArray_math.h"
#include "identityOperator.h"

namespace Gadgetron{
  
  /** \class hoIdentityOperator
      \brief Instantiation of the identity operator on the cpu.
      
      The class hoIdentityOperator is a convienience wrapper for the device independent identityOperator.
      hoIdentityOperator instantiates the identityOperator for type hoNDArray<T>.
  */
  template <class T> class hoIdentityOperator : public identityOperator< hoNDArray<T> >
  {
  public:    
    hoIdentityOperator() : identityOperator< hoNDArray<T> >() {}
    virtual ~hoIdentityOperator() {}
  }; 
}
