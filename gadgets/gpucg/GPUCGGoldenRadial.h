#ifndef GPUCGGOLDENRADIALGADGET_H
#define GPUCGGOLDENRADIALGADGET_H

#include "GPUCGGadget.h"

class GPUCGGoldenRadialGadget : public GPUCGGadget
{

 public:
  
 GPUCGGoldenRadialGadget(bool pass_on_data = false, int slice_no = 0) 
   : GPUCGGadget(pass_on_data, slice_no)
    {}


 protected:
  virtual int calculate_trajectory();
  virtual int calculate_density_compensation();
};


#endif //GPUCGGOLDENRADIALGADGET_H
