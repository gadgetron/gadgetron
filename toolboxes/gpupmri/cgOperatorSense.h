#pragma once

#include "gadgetron_export.h"
#include "cuCGMatrixOperator.h"
#include "cuNDArray.h"
#include "vector_td.h"

#include <boost/smart_ptr.hpp>

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cgOperatorSense : public cuCGMatrixOperator<REAL,typename complext<REAL>::Type>
{

public:

  typedef typename complext<REAL>::Type _complext;
  
  cgOperatorSense( int device = -1 ) : cuCGMatrixOperator<REAL,_complext>(device), ncoils_(0) {}
  virtual ~cgOperatorSense() {}

  virtual int set_csm( boost::shared_ptr< cuNDArray<_complext> > csm );

  virtual int mult_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false ) = 0;
  virtual int mult_MH( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false ) = 0;
  virtual int mult_MH_M( cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false );

  virtual void set_device( int device ){
    std::cerr << "cgOperatorSense::set_device: device needs to be set in the constructor (for now)." << std::endl; }

  virtual int set_device(){ return cuCGMatrixOperator<REAL,_complext>::set_device(); }

protected:

  int mult_csm( cuNDArray<_complext>* in, cuNDArray<_complext>* out );
  int mult_csm_conj_sum( cuNDArray<_complext>* in, cuNDArray<_complext>* out) ;

  boost::shared_ptr< cuNDArray<_complext> > csm_;
  unsigned int ncoils_;
  std::vector<unsigned int> dimensionsI_;
  std::vector<unsigned int> dimensionsK_;
};
