/** \file cuNonCartesianSenseOperator.h
    \brief Non-Cartesian Sense operator, GPU based.
*/

#pragma once

#include "cuSenseOperator.h"
#include "cuNFFT.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> class EXPORTGPUPMRI cuNonCartesianTSenseOperator : public cuSenseOperator<REAL,D>
  {
  
  public:
  
    typedef typename uint64d<D>::Type _uint64d;
    typedef typename reald<REAL,D>::Type _reald;

    cuNonCartesianTSenseOperator(ConvolutionType conv = ConvolutionType::STANDARD);
    virtual ~cuNonCartesianTSenseOperator() {}
    
    inline std::vector<boost::shared_ptr< cuNFFT_plan<REAL, D> >> get_plan() { return plan_; }
    inline std::vector<boost::shared_ptr< cuNDArray<REAL> >> get_dcw() { return dcw_; }
    inline bool is_preprocessed() { return is_preprocessed_; } 

    virtual void mult_M( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate = false );
    virtual void mult_MH( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate = false );

    virtual void setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W );
    virtual void preprocess(std::vector<cuNDArray<_reald>> &trajectory);
    virtual void set_dcw(std::vector<boost::shared_ptr<cuNDArray<REAL>>> dcw);
    virtual void set_shots_per_time(std::vector<size_t> shots_per_time);


  
  protected:
    std::vector<boost::shared_ptr< cuNFFT_plan<REAL, D> > >plan_;
    std::vector<boost::shared_ptr< cuNDArray<REAL> > > dcw_;
    std::vector<size_t> shots_per_time_;
    ConvolutionType convolutionType;
    bool is_preprocessed_;
  };
  
}
