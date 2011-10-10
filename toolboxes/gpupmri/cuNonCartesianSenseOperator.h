#pragma once

#include "cuSenseOperator.h"
#include "NFFT.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cuNonCartesianSenseOperator : public cuSenseOperator<REAL,D>
{

 public:
  
  cuNonCartesianSenseOperator( int device = -1 ) : cuSenseOperator<REAL,D>(device) {}
  virtual ~cuNonCartesianSenseOperator() {}

  typedef typename cuSenseOperator<REAL,D>::_complext _complext;
  typedef typename uintd<D>::Type _uintd;
  typedef typename reald<REAL,D>::Type _reald;
  
  virtual int mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );
  virtual int mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );

  virtual int setup( _uintd matrix_size, _uintd matrix_size_os, REAL W );
  virtual int preprocess( cuNDArray<_reald> *trajectory );
  virtual int set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw );

 protected:

  NFFT_plan<REAL, D> plan_;
  boost::shared_ptr< cuNDArray<REAL> > dcw_;
};
