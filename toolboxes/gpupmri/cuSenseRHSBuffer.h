#pragma once

#include "gadgetron_export.h"
#include "cuSenseOperator.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cuSenseRHSBuffer : public cuSenseOperator<REAL,D>
{
 public:
  
  typedef typename cuSenseOperator<REAL,D>::_complext _complext;
  
  cuSenseRHSBuffer( int device = -1 ) : cuSenseOperator<REAL,D>(device) {}
  virtual ~cuSenseRHSBuffer() {}

  virtual int mult_M( cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate = false );
  virtual int mult_MH( cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate = false );
};
