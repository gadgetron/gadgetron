#ifndef GPUCGGOLDENRADIALGADGET_H
#define GPUCGGOLDENRADIALGADGET_H
#pragma once

#include "GPUCGGadget.h"

class EXPORTGADGETSCGSENSE GPUCGGoldenRadialGadget : public GPUCGGadget
{

 public:
  GADGET_DECLARE(GPUCGGoldenRadialGadget);

 protected:
  virtual boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory();
  virtual boost::shared_ptr< cuNDArray<float> > calculate_density_compensation();
};


#endif //GPUCGGOLDENRADIALGADGET_H
