#pragma once

#include "hoCuNDArray.h"
#include "hoNDArray_utils.h"
#include "complext.h"

namespace Gadgetron{

  /**
   * @brief Calculates the elementwise absolute value of the array
   * @param[in] data Input data
   * @return A new array containing the elementwise absolute value of data
   */
  template<class T>
  boost::shared_ptr<hoCuNDArray<typename realType<T>::type> > abs(hoCuNDArray<T> *data){
    return boost::static_pointer_cast<hoCuNDArray<typename realType<T>::type> >(abs(static_cast<hoNDArray<T>* >(data)));
  }
}
