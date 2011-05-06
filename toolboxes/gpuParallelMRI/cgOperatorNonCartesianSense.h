#pragma once

#include "cgOperatorSense.h"
#include "NFFT.h"

template<class REAL, unsigned int D>
class cgOperatorNonCartesianSense : public cgOperatorSense<REAL,D>
{
 public:
  
  cgOperatorNonCartesianSense() : cgOperatorSense<REAL,D>(), trajectory_(0), weights_(0) {}

  typedef typename cgOperatorSense<REAL,D>::_complext _complext;
  typedef typename uintd<D>::Type _uintd;
  typedef typename reald<REAL,D>::Type _reald;
  
  virtual int mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );
  virtual int mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );

  virtual int setup( _uintd matrix_size, _uintd matrix_size_os, REAL W );
  virtual int set_trajectory( cuNDArray<_reald> *trajectory );
  virtual int set_weights( cuNDArray<REAL> *w );

 protected:
  cuNDArray<_reald> *trajectory_;
  cuNDArray<float> *weights_;
  NFFT_plan<REAL, D> plan_;
};
