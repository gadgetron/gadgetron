#pragma once

#include "gpuRadialSensePrepBase.h"

namespace Gadgetron{
  
  class EXPORTGADGETS_GPURADIAL gpuFixedRadialSensePrepGadget : public gpuRadialSensePrepBase
  {
  public:
    
    GADGET_DECLARE(gpuFixedRadialSensePrepGadget);
    
    gpuFixedRadialSensePrepGadget() : gpuRadialSensePrepBase() {};
    virtual ~gpuFixedRadialSensePrepGadget() {}

  protected:

    virtual boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory_for_buffer(long profile_offset);
    virtual boost::shared_ptr< hoNDArray<floatd2> > calculate_trajectory_for_reconstruction(long profile_offset);
    virtual boost::shared_ptr< cuNDArray<float> >   calculate_density_compensation_for_buffer();
    virtual boost::shared_ptr< hoNDArray<float> >   calculate_density_compensation_for_reconstruction();    
  };
}
