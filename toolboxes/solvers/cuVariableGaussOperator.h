#pragma once

#include "cuNDArray.h"
#include "linearOperator.h"
#include "cuLinearOperator_macros.h"
#include "vector_td_utilities.h"
#include "solvers_export.h"

#include <boost/smart_ptr.hpp>
#include <vector>

template <class REAL, class T, unsigned int D> class EXPORTSOLVERS cuVariableGaussOperator 
	: public linearOperator<REAL, cuNDArray<T> >
{

 public:
  cuVariableGaussOperator() : linearOperator<REAL, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuVariableGaussOperator() {}

  void set_sigma(cuNDArray<REAL> * sigma);
  
  virtual int mult_MH_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false );
  virtual int mult_MH( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false );
  virtual int mult_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false );

  virtual boost::shared_ptr< linearOperator< REAL, cuNDArray<T> > > clone(){
    return linearOperator< REAL, cuNDArray<T> >::clone(this);
  }
  
 protected:
  cuNDArray<REAL>* _sigma;
  boost::shared_ptr<cuNDArray<REAL> > _norm;

  DECLARE_LINEAR_OPERATOR_DEVICE_SUPPORT(cuVariableGaussOperator)
};
