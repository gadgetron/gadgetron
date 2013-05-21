#pragma once

#include "cuNDArray.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "resampleOperator.h"
#include "gpureg_export.h"

#include <thrust/device_vector.h>

namespace Gadgetron{

  template <class T, unsigned int D>
  class EXPORTGPUREG cuResampleOperator : public resampleOperator< cuNDArray<typename realType<T>::Type>, cuNDArray<T> >
  {    
  public:

    typedef typename realType<T>::Type REAL;
    
    cuResampleOperator() : resampleOperator< cuNDArray<REAL>, cuNDArray<T> >() {}
    virtual ~cuResampleOperator() {}
  
    virtual void reset()
    {
      lower_bounds_ = thrust::device_vector<unsigned int>();
      upper_bounds_ = thrust::device_vector<unsigned int>();
      indices_ = thrust::device_vector<unsigned int>();
      weights_ = thrust::device_vector<REAL>();
      resampleOperator< cuNDArray<typename realType<T>::Type>, cuNDArray<T> >::reset();
    }
    
    virtual void mult_MH_preprocess();
  
  protected:
    virtual unsigned int get_num_neighbors() = 0;
    virtual void write_sort_arrays( thrust::device_vector<unsigned int> &sort_keys ) = 0;
  
  protected:
    virtual void setup_grid( dim3 *blockDim, dim3* gridDim, unsigned int number_of_elements, unsigned int num_batches = 1 );
  
  protected:
    thrust::device_vector<unsigned int> lower_bounds_;
    thrust::device_vector<unsigned int> upper_bounds_;
    thrust::device_vector<unsigned int> indices_;
    thrust::device_vector<REAL> weights_;
  };
}
