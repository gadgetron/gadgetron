#pragma once

#include "senseOperator.h"
#include "cuLinearOperator_macros.h"
#include "gpupmri_export.h"
#include "cuNDArray.h"
#include "vector_td.h"


template<class REAL, unsigned int D> class EXPORTGPUPMRI cuSenseOperator 
	: public senseOperator<REAL, D, cuNDArray< complext<REAL> > >
{

public:

  cuSenseOperator() : senseOperator<REAL,D,cuNDArray< _complext> >() { set_device(-1); }
  virtual ~cuSenseOperator() {}

  typedef complext<REAL> _complext;

  virtual void mult_M( cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate = false ) = 0;
  virtual void mult_MH( cuNDArray<_complext> *in, cuNDArray<_complext> *out, bool accumulate = false ) = 0;
    
  virtual void mult_csm( cuNDArray<_complext> *in, cuNDArray< _complext> *out );
  virtual void mult_csm_conj_sum( cuNDArray<_complext> *in, cuNDArray<_complext> *out );
  
  DECLARE_LINEAR_OPERATOR_DEVICE_SUPPORT(cuSenseOperator)
};
