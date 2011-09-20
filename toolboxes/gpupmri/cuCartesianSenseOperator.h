#pragma once

#include "cuSenseOperator.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cuCartesianSenseOperator : public cuSenseOperator<REAL,D>
{
 public:

  cuCartesianSenseOperator( int device = -1 ) : cuSenseOperator<REAL,D>(device) {}
  virtual ~cuCartesianSenseOperator() {}

  typedef typename cuSenseOperator<REAL,D>::_complext _complext;

  virtual int mult_M( cuNDArray< _complext> *in, cuNDArray< _complext> *out, bool accumulate = false);
  virtual int mult_MH( cuNDArray< _complext> *in, cuNDArray< _complext> *out, bool accumulate = false);

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
