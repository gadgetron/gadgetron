/** \file cuNonCartesianSenseOperator.h
    \brief Non-Cartesian Sense operator, GPU based.
*/

#pragma once

#include "cuSenseOperator.h"
#include "cuNFFT.h"

namespace Gadgetron{

  template<class REAL, unsigned int D, bool ATOMICS = false> class EXPORTGPUPMRI cuNonCartesianSenseOperator : public cuSenseOperator<REAL,D>
  {
  
  public:
  
    typedef typename uint64d<D>::Type _uint64d;
    typedef typename reald<REAL,D>::Type _reald;

    cuNonCartesianSenseOperator() : cuSenseOperator<REAL,D>() { 
      plan_ = boost::shared_ptr< cuNFFT_plan<REAL, D, ATOMICS> >( new cuNFFT_plan<REAL, D, ATOMICS>() );
      is_preprocessed_ = false;
    }
    
    virtual ~cuNonCartesianSenseOperator() {}
    
    inline boost::shared_ptr< cuNFFT_plan<REAL, D, ATOMICS> > get_plan() { return plan_; }
    inline boost::shared_ptr< cuNDArray<REAL> > get_dcw() { return dcw_; }
    inline bool is_preprocessed() { return is_preprocessed_; } 

    virtual void mult_M( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate = false );
    virtual void mult_MH( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate = false );

    virtual void setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W );
    virtual void preprocess( cuNDArray<_reald> *trajectory );
    virtual void set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw );


  
  protected:
    boost::shared_ptr< cuNFFT_plan<REAL, D, ATOMICS> > plan_;
    boost::shared_ptr< cuNDArray<REAL> > dcw_;
    bool is_preprocessed_;
  };
  
  //Atomics can't be used with doubles
  template<unsigned int D> class EXPORTGPUPMRI cuNonCartesianSenseOperator<double,D,true>{};
}
