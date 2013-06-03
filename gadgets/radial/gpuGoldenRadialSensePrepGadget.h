#pragma once

#include "gpuRadialSensePrepBase.h"

namespace Gadgetron{
  
  class EXPORTGADGETS_GPURADIAL gpuGoldenRadialSensePrepGadget : public gpuRadialSensePrepBase
  {
  public:
    
    GADGET_DECLARE(gpuGoldenRadialSensePrepGadget);
    
    gpuGoldenRadialSensePrepGadget();
    virtual ~gpuGoldenRadialSensePrepGadget() {}

  protected:

    virtual int process_config(ACE_Message_Block* mb);

    virtual boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory_for_buffer(long profile_offset);
    virtual boost::shared_ptr< hoNDArray<floatd2> > calculate_trajectory_for_reconstruction(long profile_offset);
    virtual boost::shared_ptr< cuNDArray<float> >   calculate_density_compensation_for_buffer();
    virtual boost::shared_ptr< hoNDArray<float> >   calculate_density_compensation_for_reconstruction();  
  };
}
