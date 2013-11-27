/** \file cuUpsampleOperator.h
    \brief Instantiation of the upsampling operator on the gpu.
    
    The file cuUpsampleOperator.h is a convienience wrapper for the device independent upsampleOperator class.
    The class cuUpsampleOperator instantiates the upsampleOperator for the cuNDArray
    and the header furthermore includes additional neccessary header files.
*/

#pragma once

#include "cuNDArray_utils.h"
#include "upsampleOperator.h"

namespace Gadgetron{
  
  /** \class cuUpsampleOperator
      \brief Instantiation of the upsample operator on the gpu.
      
      The class cuUpsampleOperator is a convienience wrapper for the device independent upsampleOperator.
      cuUpsampleOperator instantiates the upsampleOperator for type cuNDArray<T>.
  */
  template <class T, unsigned int D> class cuUpsampleOperator : public upsampleOperator<cuNDArray<T>, D>
  {
  public:    
    cuUpsampleOperator() : upsampleOperator<cuNDArray<T>,D>() {}
    virtual ~cuUpsampleOperator() {}
  }; 
}
