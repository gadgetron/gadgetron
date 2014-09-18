#pragma once

#include "hoNDArray.h"
#include "cpucore_math_export.h"

namespace Gadgetron{


	/***
	 * Finds the maximum element of the array
	 */
  template<class REAL> EXPORTCPUCOREMATH REAL max(hoNDArray<REAL>* data);

  /***
   * Finds the minimum element of the array
   */
  template<class REAL> EXPORTCPUCOREMATH REAL min(hoNDArray<REAL>* data);

  /***
   * Finds the mean of the array
   */
  template<class T> EXPORTCPUCOREMATH T mean(hoNDArray<T>* data);

  /***
   * Calculates the sum of the array
   */
  template<class T> EXPORTCPUCOREMATH T sum(hoNDArray<T>* data);

  /***
   * Calculates the std of the array
   */
  template<class T> EXPORTCPUCOREMATH T stddev(hoNDArray<T>* data);
}
