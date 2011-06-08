#ifndef GPUCGGOLDENRADIALGADGET_H
#define GPUCGGOLDENRADIALGADGET_H
#pragma once

#include "GPUCGGadget.h"

class EXPORTGADGETSCGSENSE GPUCGGoldenRadialGadget : public GPUCGGadget
{

 public:
  GADGET_DECLARE(GPUCGGoldenRadialGadget);

 protected:
  virtual int calculate_trajectory();
  virtual int calculate_density_compensation();
};


#endif //GPUCGGOLDENRADIALGADGET_H
