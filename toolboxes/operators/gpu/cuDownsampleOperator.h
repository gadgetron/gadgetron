/** \file cuDownsampleOperator.h
    \brief Instantiation of the downsampling operator on the gpu.
    
    The file cuDownsampleOperator.h is a convienience wrapper for the device independent downsampleOperator class.
    The class cuDownsampleOperator instantiates the downsampleOperator for the cuNDArray
    and the header furthermore includes additional neccessary header files.
*/

#pragma once

#include "cuNDArray_utils.h"
#include "downsampleOperator.h"

namespace Gadgetron{
  
  /** \class cuDownsampleOperator
      \brief Instantiation of the downsample operator on the gpu.
      
      The class cuDownsampleOperator is a convienience wrapper for the device independent downsampleOperator.
      cuDownsampleOperator instantiates the downsampleOperator for type cuNDArray<T>.
  */
  template <class T, unsigned int D> class cuDownsampleOperator : public downsampleOperator<cuNDArray<T>,D>
  {
  public:    
    cuDownsampleOperator() : downsampleOperator<cuNDArray<T>,D>() {}
    virtual ~cuDownsampleOperator() {}
  }; 
}
