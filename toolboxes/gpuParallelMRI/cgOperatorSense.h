#pragma once

#include "cuCGMatrixOperator.h"
#include "cuNDArray.h"
#include "vector_td.h"

template<class REAL, unsigned int D>
class cgOperatorSense : public cuCGMatrixOperator< typename complext<REAL>::Type >
{
public:

  cgOperatorSense() : csm_(0x0), ncoils_(0) {}

  typedef typename complext<REAL>::Type _complext;
  
  virtual int set_csm( cuNDArray<_complext>* csm );

  virtual int mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false ) = 0;
  virtual int mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false ) = 0;
  virtual int mult_MH_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );

protected:
  cuNDArray<_complext>* csm_;
  unsigned int ncoils_;
  std::vector<unsigned int> dimensionsI_;
  std::vector<unsigned int> dimensionsK_;
  
  int mult_csm( cuNDArray<_complext>* in, cuNDArray<_complext>* out );
  int mult_csm_conj_sum( cuNDArray<_complext>* in, cuNDArray<_complext>* out) ;
};
