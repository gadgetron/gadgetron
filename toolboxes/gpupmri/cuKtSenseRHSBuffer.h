#pragma once

#include "cuSenseRHSBuffer.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cuKtSenseRHSBuffer : public cuSenseRHSBuffer<REAL,D>
{
 public:
  
  typedef typename cuSenseOperator<REAL,D>::_complext _complext;
  
  cuKtSenseRHSBuffer( int device = -1 ) : cuSenseRHSBuffer<REAL,D>(device) {}
  virtual ~cuKtSenseRHSBuffer() {}

  virtual int mult_MH( cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate = false );
};
