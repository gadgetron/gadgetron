#pragma once

#include "cgOperatorSenseRHSBuffer.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cgOperatorKtSenseRHSBuffer : public cgOperatorSenseRHSBuffer<REAL,D>
{
 public:
  
  typedef typename cgOperatorSense<REAL,D>::_complext _complext;
  
  cgOperatorKtSenseRHSBuffer( int device = -1 ) : cgOperatorSenseRHSBuffer<REAL,D>(device) {}
  virtual ~cgOperatorKtSenseRHSBuffer() {}

  virtual int mult_MH( cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate = false );
};
