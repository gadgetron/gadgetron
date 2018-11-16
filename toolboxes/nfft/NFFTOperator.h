#pragma once

#include "nfft_export.h"
#include "linearOperator.h"
#include "NFFT.h"
#include "vector_td.h"

namespace Gadgetron{

  template<template<class> class ARRAY, class REAL, unsigned int D> 
  class EXPORTNFFT NFFTOperator : public virtual linearOperator<ARRAY< complext<REAL> > >
  {  
  public:


    NFFTOperator();
  
    virtual ~NFFTOperator() {}
  
    virtual void set_dcw( boost::shared_ptr< ARRAY<REAL> > dcw ) { dcw_ = dcw; }
    inline boost::shared_ptr< ARRAY<REAL> > get_dcw() { return dcw_; }

    inline boost::shared_ptr<NFFT_plan<ARRAY,REAL,D>> get_plan() { return plan_; }
  
    virtual void setup( typename uint64d<D>::Type matrix_size, typename uint64d<D>::Type matrix_size_os, REAL W );
    virtual void preprocess(const ARRAY<typename reald<REAL,D>::Type>& trajectory );

    virtual void mult_M( ARRAY< complext<REAL> > *in, ARRAY< complext<REAL> > *out, bool accumulate = false );
    virtual void mult_MH( ARRAY< complext<REAL> > *in, ARRAY< complext<REAL> > *out, bool accumulate = false );
    virtual void mult_MH_M( ARRAY< complext<REAL> > *in, ARRAY< complext<REAL> > *out, bool accumulate = false );


  protected:
    boost::shared_ptr<NFFT_plan<ARRAY,REAL,D>> plan_;
    boost::shared_ptr< ARRAY<REAL> > dcw_;
  };
}
