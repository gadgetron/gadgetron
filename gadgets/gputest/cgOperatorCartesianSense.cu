#include "cgOperatorCartesianSense.h"

int cgOperatorCartesianSense::mult_M(cuNDArray<float2>* in, 
				     cuNDArray<float2>* out, 
				     bool accumulate)
{
  if (!(in->dimensions_equal(dimensions_)) ||
      !(out->dimensions_equal(dimensions_out_)) ) {
    std::cerr << "cgOperatorCartesianSense::mult_M dimensions mismatch" << std::endl;
    return -1;
  }


  return 0;
}

int cgOperatorCartesianSense::mult_MH(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate)
{

  if (!(out->dimensions_equal(dimensions_)) ||
      !(in->dimensions_equal(dimensions_out_)) ) {
    std::cerr << "cgOperatorCartesianSense::mult_MH dimensions mismatch" << std::endl;
    return -1;
  }


  return 0;
}

int cgOperatorCartesianSense::mult_MH_M(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate)
{
  /*
  cuNDArray<float2> tmp;
  if (!tmp->create(dimensions_out_)) {
    std::cerr << "cgOperatorCartesianSense::mult_MH_M: Unable to create temporary storage" << std::endl;
    return -1;
  }

  if (mult_M(in, &tmp, false) < 0) {
    std::cerr << "cgOperatorCartesianSense::mult_MH_M: Unable to perform mult_M" << std::endl;
    return -2;
  }

  if (mult_M(&tmp, out, accumulate) < 0) {
    std::cerr << "cgOperatorCartesianSense::mult_MH_M: Unable to perform mult_M" << std::endl;
    return -2;
  }
  */

  return 0;
}
