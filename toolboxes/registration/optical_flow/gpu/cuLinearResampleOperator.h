#pragma once

#include "cuResampleOperator.h"

template <class REAL, class T, unsigned int D>
class EXPORTGPUREG cuLinearResampleOperator : public cuResampleOperator<REAL,T,D>
{
  
 public:
  
  cuLinearResampleOperator() : cuResampleOperator<REAL,T,D>() {}
  virtual ~cuLinearResampleOperator() {}
  
  virtual int mult_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false);
  virtual int mult_MH( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false);
  virtual int mult_MH_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false);
  
  virtual boost::shared_ptr< linearOperator< REAL, cuNDArray<T> > > clone() {
    return linearOperator< REAL, cuNDArray<T> >::clone(this);
  }
  
 protected:
  virtual unsigned int get_num_neighbors();
  virtual bool write_sort_arrays( thrust::device_vector<unsigned int> &sort_keys );
};
