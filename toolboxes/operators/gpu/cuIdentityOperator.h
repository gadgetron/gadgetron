/** \file cuIdentityOperator.h
    \brief Instantiation of the identity operator on the gpu.
    
    The file cuIdentityOperator.h is a convienience wrapper for the device independent identityOperator class.
    The class cuIdentityOperator instantiates the identityOperator for the cuNDArray
    and the header furthermore includes additional neccessary header files.
*/

#pragma once

#include "cuNDArray_math.h"
#include "identityOperator.h"

namespace Gadgetron{
  
  /** \class cuIdentityOperator
      \brief Instantiation of the identity operator on the gpu.
      
      The class cuIdentityOperator is a convienience wrapper for the device independent identityOperator.
      cuIdentityOperator instantiates the identityOperator for type cuNDArray<T>.
  */
  template <class T> class cuIdentityOperator : public identityOperator< cuNDArray<T> >
  {
  public:    
    cuIdentityOperator() : identityOperator< cuNDArray<T> >() {}
    virtual ~cuIdentityOperator() {}
  }; 
}
