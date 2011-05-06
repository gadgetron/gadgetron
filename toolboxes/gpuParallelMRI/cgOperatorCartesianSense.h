#pragma once

#include "cgOperatorSense.h"

template<class REAL, unsigned int D>
class cgOperatorCartesianSense : public cgOperatorSense<REAL,D>
{
 public:

  cgOperatorCartesianSense() : cgOperatorSense<REAL,D>(), idx_(0) {}

  typedef typename cgOperatorSense<REAL,D>::_complext _complext;

  virtual int mult_M(cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false);
  virtual int mult_MH(cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false);
  //virtual int mult_MH_M(cuNDArray<_complext>* in, cuNDArray<_complext>* out, bool accumulate = false);

  virtual int set_csm(cuNDArray<_complext>* csm) {
    if (csm != 0) {
      this->csm_ = csm;
      this->ncoils_ = csm->get_size(csm->get_number_of_dimensions()-1);
      this->dimensions_ = csm->get_dimensions();
      this->dimensions_.pop_back();
    }
    return 0;
  }

  virtual int set_sampling_indices(cuNDArray<unsigned int>* idx) {
    if (idx) {
      idx_ = idx;
      this->nsamples_ = idx_->get_number_of_elements();
      this->dimensions_out_.clear();
      this->dimensions_out_.push_back(this->nsamples_);
      this->dimensions_out_.push_back(this->ncoils_);
    }
    return 0;
  }

 protected:
  cuNDArray<unsigned int>* idx_;
};
