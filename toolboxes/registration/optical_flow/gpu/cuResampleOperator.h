#pragma once

#include "resampleOperator.h"
#include "cuLinearOperator_macros.h"
#include "cuNDArray.h"
#include "gpureg_export.h"

#include <thrust/device_vector.h>

template <class REAL, class T, unsigned int D>
class EXPORTGPUREG cuResampleOperator : public resampleOperator< REAL, cuNDArray<REAL>, cuNDArray<T> >
{
  
 public:
  
  cuResampleOperator() : resampleOperator< REAL, cuNDArray<REAL>, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuResampleOperator() {}
  
  virtual void reset(){

    lower_bounds_ = thrust::device_vector<unsigned int>();
    upper_bounds_ = thrust::device_vector<unsigned int>();
    indices_ = thrust::device_vector<unsigned int>();
    weights_ = thrust::device_vector<REAL>();

    resampleOperator< REAL, cuNDArray<REAL>, cuNDArray<T> >::reset();
  }

  virtual bool mult_MH_preprocess();
  
 protected:
  virtual unsigned int get_num_neighbors() = 0;
  virtual bool write_sort_arrays( thrust::device_vector<unsigned int> &sort_keys ) = 0;
  
 protected:
  virtual bool setup_grid( dim3 *blockDim, dim3* gridDim, unsigned int number_of_elements, unsigned int num_batches = 1 );
  
 protected:
  thrust::device_vector<unsigned int> lower_bounds_;
  thrust::device_vector<unsigned int> upper_bounds_;
  thrust::device_vector<unsigned int> indices_;
  thrust::device_vector<REAL> weights_;
  
  DECLARE_LINEAR_OPERATOR_DEVICE_SUPPORT(cuResampleOperator)
};
