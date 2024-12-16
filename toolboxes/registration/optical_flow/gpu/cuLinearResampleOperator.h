#pragma once

#include "cuResampleOperator.h"

namespace Gadgetron{

  template <class T, unsigned int D>
  class cuLinearResampleOperator : public cuResampleOperator<T,D>
  {
  public:

    cuLinearResampleOperator() : cuResampleOperator<T,D>() {}
    virtual ~cuLinearResampleOperator() {}

    virtual void mult_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false);
    virtual void mult_MH( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false);


  protected:
    virtual unsigned int get_num_neighbors();
    virtual void write_sort_arrays( thrust::device_vector<unsigned int> &sort_keys );
  };
}
