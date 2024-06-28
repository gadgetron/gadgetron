/** \file cuCartesianSenseOperator.h
    \brief Cartesian Sense operator, GPU based.
*/

#pragma once

#include "cuSenseOperator.h"

namespace Gadgetron{
  
  template<class REAL, unsigned int D> class EXPORTGPUPMRI cuCartesianSenseOperator : public cuSenseOperator<REAL,D>
  {
  public:
    
    cuCartesianSenseOperator() : cuSenseOperator<REAL,D>() {}
    virtual ~cuCartesianSenseOperator() {}
    
    virtual void mult_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false);
    virtual void mult_MH( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false);
    
    virtual void set_sampling_indices( boost::shared_ptr< cuNDArray<unsigned int> > idx) 
    {
      if (idx.get()) {
        idx_ = idx;
        std::vector<size_t> tmp_dims;
        tmp_dims.push_back(idx_->get_number_of_elements());
        tmp_dims.push_back(this->ncoils_);
        this->set_codomain_dimensions(tmp_dims);
      }
    }
    
  protected:
    boost::shared_ptr< cuNDArray<unsigned int> > idx_;
  };
}
