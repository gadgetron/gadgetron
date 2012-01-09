#ifndef GPUCGFIXEDRADIALGADGET_H
#define GPUCGFIXEDRADIALGADGET_H
#pragma once

#include "GPUCGGadget.h"

class EXPORTGADGETSCGSENSE GPUCGFixedRadialGadget : public GPUCGGadget
{

 public:
  GADGET_DECLARE(GPUCGFixedRadialGadget);

 protected:
  virtual boost::shared_ptr< cuNDArray<floatd2::Type> > calculate_trajectory();
  virtual boost::shared_ptr< cuNDArray<float> > calculate_density_compensation();
};


#endif //GPUCGGOLDENRADIALGADGET_H
