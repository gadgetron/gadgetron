#pragma once

#include "gadgetron_export.h"
#include "cgOperatorSense.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cgOperatorSenseRHSBuffer : public cgOperatorSense<REAL,D>
{
 public:
  
  typedef typename cgOperatorSense<REAL,D>::_complext _complext;
  
  cgOperatorSenseRHSBuffer( int device = -1 ) : cgOperatorSense<REAL,D>(device) {}
  virtual ~cgOperatorSenseRHSBuffer() {}

  virtual int mult_M( cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate = false );
  virtual int mult_MH( cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate = false );
};
