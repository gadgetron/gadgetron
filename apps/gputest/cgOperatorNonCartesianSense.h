#ifndef CGOPERATORNONCARTESIANSENSE_H
#define CGOPERATORNONCARTESIANSENSE_H

#include "cgOperatorSense.h"
#include "NFFT.h"

class cgOperatorNonCartesianSense : public cgOperatorSense
{

 public:
  
 cgOperatorNonCartesianSense()
   : cgOperatorSense()
    , trajectory_(0)
    , weights_(0)
    { }

  virtual int mult_M(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate = false);
  virtual int mult_MH(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate = false);
  virtual int set_trajectories(cuNDArray<floatd2>* trajectory);

  virtual int set_weights(cuNDArray<float>* w) {
     if (!trajectory_) {
      std::cerr << "cgOperatorNonCartesianSense::set_weights : Error setting weights, trajectory not set" << std::endl;
      return -1;
    }
    if (w->get_number_of_elements() != trajectory_->get_number_of_elements()) {
      std::cerr << "cgOperatorNonCartesianSense::set_weights : weights don't match trajectory" << std::endl;
       return -1;
    }

    weights_ = w;

    return 0;
  }

 protected:
  cuNDArray<floatd2>* trajectory_;
  cuNDArray<float>* weights_;
  NFFT_plan<float, 2> plan_;
  
};

#endif //CGOPERATORNONCARTESIANSENSE_H
