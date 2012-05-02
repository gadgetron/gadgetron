#pragma once

#include "cuNonCartesianSenseOperator.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cuNonCartesianKtSenseOperator : public cuNonCartesianSenseOperator<REAL,D>
{

 public:
  
  typedef typename cuSenseOperator<REAL,D>::_complext _complext;
  typedef typename uintd<D>::Type _uintd;
  typedef typename reald<REAL,D>::Type _reald;

  cuNonCartesianKtSenseOperator() : cuNonCartesianSenseOperator<REAL,D>() {}
  virtual ~cuNonCartesianKtSenseOperator() {}
  
  virtual int mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );
  virtual int mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );

  virtual boost::shared_ptr< linearOperator< REAL, cuNDArray< complext<REAL>  > > > clone(){
    return linearOperator< REAL, cuNDArray<complext<REAL> > >::clone(this);
  }  
};
