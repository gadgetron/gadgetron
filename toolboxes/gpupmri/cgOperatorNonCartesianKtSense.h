#pragma once

#include "cgOperatorNonCartesianSense.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cgOperatorNonCartesianKtSense : public cgOperatorNonCartesianSense<REAL,D>
{

 public:
  
  cgOperatorNonCartesianKtSense( int device = -1 ) : cgOperatorNonCartesianSense<REAL,D>(device) {}
  virtual ~cgOperatorNonCartesianKtSense() {}

  typedef typename cgOperatorSense<REAL,D>::_complext _complext;
  typedef typename cgOperatorNonCartesianSense<REAL,D>::_uintd _uintd;
  typedef typename cgOperatorNonCartesianSense<REAL,D>::_reald _reald;
  
  virtual int mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );
  virtual int mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );

};
