#pragma once

#include "cuNonCartesianSenseOperator.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cuNonCartesianKtSenseOperator : public cuNonCartesianSenseOperator<REAL,D>
{

 public:
  
  cuNonCartesianKtSenseOperator( int device = -1 ) : cuNonCartesianSenseOperator<REAL,D>(device) {}
  virtual ~cuNonCartesianKtSenseOperator() {}

  typedef typename cuSenseOperator<REAL,D>::_complext _complext;
  typedef typename cuNonCartesianSenseOperator<REAL,D>::_uintd _uintd;
  typedef typename cuNonCartesianSenseOperator<REAL,D>::_reald _reald;
  
  virtual int mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );
  virtual int mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );

};
