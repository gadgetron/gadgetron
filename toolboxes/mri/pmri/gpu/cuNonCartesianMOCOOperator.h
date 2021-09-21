/** \file cuNonCartesianSenseOperator.h
    \brief Non-Cartesian Sense operator, GPU based.
*/

#pragma once

#include "cuSenseOperator.h"
#include "cuNFFT.h"
#include "cuCKOpticalFlowSolver.h"
#include "cuLinearResampleOperator.h"
#include "hoNDArray.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> class EXPORTGPUPMRI cuNonCartesianMOCOOperator : public cuSenseOperator<REAL,D>
  {
  
  public:
  
    typedef typename uint64d<D>::Type _uint64d;
    typedef typename reald<REAL,D>::Type _reald;

    cuNonCartesianMOCOOperator(ConvolutionType conv = ConvolutionType::STANDARD);
    virtual ~cuNonCartesianMOCOOperator() {}
    
    inline std::vector<boost::shared_ptr< cuNFFT_plan<REAL, D> >> get_plan() { return plan_; }
    inline std::vector< cuNDArray<REAL> > get_dcw() { return dcw_; }
    inline bool is_preprocessed() { return is_preprocessed_; } 

    void mult_M( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate = false );
    void mult_MH( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate = false );

    virtual void setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W );
    virtual void preprocess(std::vector<cuNDArray<_reald>> &trajectory);
    virtual void set_dcw( std::vector<cuNDArray<REAL> > dcw);
    virtual void set_shots_per_time(std::vector<size_t> shots_per_time);
    
    void set_forward_deformation (std::vector<cuNDArray<REAL>> forward_deformation);
    void set_backward_deformation(std::vector<cuNDArray<REAL>> backward_deformation);
    
    private:
    void applyDeformation(cuNDArray< complext<REAL> > *moving_image, cuNDArray<REAL>  transformation);

  
  protected:
    std::vector<boost::shared_ptr< cuNFFT_plan<REAL, D> > >plan_;
    std::vector<cuNDArray<REAL> > dcw_;
    std::vector<cuNDArray<REAL> > forward_deformation_;
    std::vector<cuNDArray<REAL> > backward_deformation_;
    // nhlbi_toolbox::motion_correction moco;
    std::vector<size_t> shots_per_time_;
    ConvolutionType convolutionType;
    bool is_preprocessed_;
  };
  
}
