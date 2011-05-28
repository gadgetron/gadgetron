#pragma once

#include "cgOperatorSense.h"

template<class REAL, unsigned int D>
class cgOperatorSenseRHSBuffer : public cgOperatorSense<REAL,D>
{
 public:
  
  typedef typename cgOperatorSense<REAL,D>::_complext _complext;
  
  cgOperatorSenseRHSBuffer( std::auto_ptr< cuNDArray<_complext> > csm ) : cgOperatorSense<REAL,D>() { set_csm(csm); }

  virtual int mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );
  virtual int mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );
};
