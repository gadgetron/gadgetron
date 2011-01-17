#ifndef GPUCGGOLDENSPIRALGADGET_H
#define GPUCGGOLDENSPIRALGADGET_H

#include "GPUCGGadget.h"

class GPUCGGoldenSpiralGadget : public GPUCGGadget
{

 protected:
  virtual int calculate_trajectory();
  virtual int calculate_density_compensation();

};

#endif //GPUCGGOLDENSPIRALGADGET_H
