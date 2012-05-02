#pragma once

#include "linearOperator.h"
#include "complext.h"
#include "NFFT.h"

template<class REAL, unsigned int D> class EXPORTGPUNFFT cuNFFTOperator 
  : public linearOperator<REAL, cuNDArray< complext<REAL> > >
{
  
 public:
  
  cuNFFTOperator() : linearOperator<REAL, cuNDArray< complext<REAL> > >() { 
    plan_ = boost::shared_ptr< NFFT_plan<REAL, D> >( new NFFT_plan<REAL, D>() );
    ready_ = false; 
  }
  
  virtual ~cuNFFTOperator() {}
  
  virtual void set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw ) { dcw_ = dcw; }
  inline boost::shared_ptr< cuNDArray<REAL> > get_dcw() { return dcw_; }

  inline boost::shared_ptr< NFFT_plan<REAL, D> > get_plan() { return plan_; }
  
  virtual int setup( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, REAL W );
  virtual int preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory );

  virtual int mult_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false );
  virtual int mult_MH( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false );
  virtual int mult_MH_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false );

  virtual boost::shared_ptr< linearOperator< REAL, cuNDArray< complext<REAL>  > > > clone(){
    return linearOperator< REAL, cuNDArray<complext<REAL> > >::clone(this);
  }
  
protected:
  boost::shared_ptr< NFFT_plan<REAL, D> > plan_;
  boost::shared_ptr< cuNDArray<REAL> > dcw_;
  std::vector<unsigned int> dimensionsK_;
  bool ready_;
};
