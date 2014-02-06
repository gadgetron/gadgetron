#pragma once

#include "gpuRadialPrepGadget.h"

namespace Gadgetron{

  class EXPORTGADGETS_RADIAL gpuRadialSensePrepGadget : public gpuRadialPrepGadget
  {
    
  public:
    GADGET_DECLARE(gpuRadialSensePrepGadget);

    gpuRadialSensePrepGadget() : gpuRadialPrepGadget() {}
    virtual ~gpuRadialSensePrepGadget() {}
    
  protected:
    
    virtual void reconfigure(unsigned int set, unsigned int slice);

    virtual boost::shared_ptr< hoNDArray<float_complext> > compute_csm( unsigned int buffer_idx );

    virtual boost::shared_ptr< hoNDArray<float_complext> > compute_reg( unsigned int set, 
                                                                        unsigned int slice, 
                                                                        bool new_frame );
    
    virtual void allocate_accumulation_buffer( unsigned int num_buffers );
    
  };

  GADGET_FACTORY_DECLARE(gpuRadialSensePrepGadget)
}
