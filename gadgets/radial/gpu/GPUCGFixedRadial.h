#ifndef GPUCGFIXEDRADIALGADGET_H
#define GPUCGFIXEDRADIALGADGET_H
#pragma once

#include "GPUCGGadget.h"
namespace Gadgetron{
class EXPORTGADGETSCGSENSE GPUCGFixedRadialGadget : public GPUCGGadget
{

 public:
  GADGET_DECLARE(GPUCGFixedRadialGadget);

  GPUCGFixedRadialGadget();

  virtual int process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2 );
  virtual int process_config( ACE_Message_Block* mb );

 protected:
  virtual boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory();
  virtual boost::shared_ptr< cuNDArray<float> > calculate_density_compensation();

  unsigned int total_projections_;
  unsigned int dynamic_acceleration_factor_;
  int previous_projection_;
};
}

#endif //GPUCGGOLDENRADIALGADGET_H
