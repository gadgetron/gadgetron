#pragma once

#include "hoNDArray.h"
#include "cpucore_math_export.h"

namespace Gadgetron{

  template<class REAL> EXPORTCPUCOREMATH REAL max(hoNDArray<REAL>* data);
  template<class REAL> EXPORTCPUCOREMATH REAL min(hoNDArray<REAL>* data);
  template<class T> EXPORTCPUCOREMATH T mean(hoNDArray<T>* data);
  template<class T> EXPORTCPUCOREMATH T sum(hoNDArray<T>* data);
}
