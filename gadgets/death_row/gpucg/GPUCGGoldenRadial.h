#ifndef GPUCGGOLDENRADIALGADGET_H
#define GPUCGGOLDENRADIALGADGET_H

#include "GPUCGGadget.h"

class GPUCGGoldenRadialGadget : public GPUCGGadget
{

 public:
  GADGET_DECLARE(GPUCGGoldenRadialGadget);

 protected:
  virtual int calculate_trajectory();
  virtual int calculate_density_compensation();
};


#endif //GPUCGGOLDENRADIALGADGET_H
