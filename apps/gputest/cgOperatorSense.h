#ifndef CGOPERATORSENSE_H
#define CGOPERATORSENSE_H

#include "cuCGMatrixOperator.h"
#include "cuNDArray.h"
#include <cuComplex.h>

class cgOperatorSense : public cuCGMatrixOperator<float2>
{

 public:
  cgOperatorSense()
    : csm_(0)
    , coils_(0)
    , samples_(0)
    { }

  virtual int mult_M(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate = false) = 0;
  virtual int mult_MH(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate = false) = 0;
  virtual int mult_MH_M(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate = false);

  virtual int set_csm(cuNDArray<float2>* csm) {
    if (csm != 0) {
      csm_ = csm;
      coils_ = csm_->get_size(csm_->get_number_of_dimensions()-1);
      dimensions_ = csm->get_dimensions();
      dimensions_.pop_back();
    }
    return 0;
  }



 protected:
  cuNDArray<float2>* csm_;
  unsigned int coils_;
  unsigned int samples_;
  std::vector<unsigned int> dimensions_;
  std::vector<unsigned int> dimensions_out_;
  
  int clear(cuNDArray<float2>* in);
  int mult_csm(cuNDArray<float2>* in, cuNDArray<float2>* out);
  int mult_csm_conj_sum(cuNDArray<float2>* in, cuNDArray<float2>* out);
};


#endif //CGOPERATOSENSE_H
