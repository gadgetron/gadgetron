#pragma once

#include "linearOperator.h"
#include "cuNFFT.h"
#include "gpunfft_export.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> class EXPORTGPUNFFT cuNFFTOperator : public linearOperator<cuNDArray< complext<REAL> > >
  {  
  public:
  
    cuNFFTOperator() : linearOperator<cuNDArray< complext<REAL> > >() {
      plan_ = boost::shared_ptr< cuNFFT_plan<REAL, D> >( new cuNFFT_plan<REAL, D>() );
      ready_ = false; 
    }
  
    virtual ~cuNFFTOperator() {}
  
    virtual void set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw ) { dcw_ = dcw; }
    inline boost::shared_ptr< cuNDArray<REAL> > get_dcw() { return dcw_; }

    inline boost::shared_ptr< cuNFFT_plan<REAL, D> > get_plan() { return plan_; }
  
    virtual void setup( typename uintd<D>::Type matrix_size, typename uintd<D>::Type matrix_size_os, REAL W );
    virtual void preprocess( cuNDArray<typename reald<REAL,D>::Type> *trajectory );

    virtual void mult_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false );
    virtual void mult_MH( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false );
    virtual void mult_MH_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false );

    virtual boost::shared_ptr< linearOperator< cuNDArray< complext<REAL>  > > > clone(){
      return linearOperator< cuNDArray<complext<REAL> > >::clone(this);
    }

  protected:
    boost::shared_ptr< cuNFFT_plan<REAL, D> > plan_;
    boost::shared_ptr< cuNDArray<REAL> > dcw_;
    std::vector<unsigned int> dimensionsK_;
    bool ready_;
  };
}
