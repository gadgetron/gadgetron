#pragma once

#include "cgOperatorSense.h"

template<class REAL, unsigned int D>
class cgOperatorCartesianSense : public cgOperatorSense<REAL,D>
{
 public:

  cgOperatorCartesianSense() : cgOperatorSense<REAL,D>() {}

  typedef typename cgOperatorSense<REAL,D>::_complext _complext;

  virtual int mult_M(cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false);
  virtual int mult_MH(cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false);
  //virtual int mult_MH_M(cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false);

  virtual int set_sampling_indices( boost::shared_ptr< cuNDArray<unsigned int> > idx) {
    if (idx.get()) {
      idx_ = idx;
      this->dimensionsK_.clear();
      this->dimensionsK_.push_back(idx_->get_number_of_elements());
      this->dimensionsK_.push_back(this->ncoils_);
    }
    return 0;
  }

 protected:
  boost::shared_ptr< cuNDArray<unsigned int> > idx_;
};
