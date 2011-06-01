#ifndef CUCGMATRIXOPERATOR_H
#define CUCGMATRIXOPERATOR_H

#pragma once
#include "gadgetron_export.h"
#include "cuNDArray.h"
#include "vector_td_utilities.h"

template <class REAL, class T> class cuCGMatrixOperator
{

 public:
  cuCGMatrixOperator() { weight_ = get_one<REAL>(); }
  virtual ~cuCGMatrixOperator() {}

  inline void set_weight( REAL weight ){ weight_ = weight; }
  inline REAL get_weight(){ return weight_; }

  virtual int mult_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  virtual int mult_MH(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  virtual int mult_MH_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  
private:
  REAL weight_;
};

#endif //CUCGMATRIXOPERATOR_H
