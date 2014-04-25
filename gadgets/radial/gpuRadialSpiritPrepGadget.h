#pragma once

#include "gpuRadialPrepGadget.h"

namespace Gadgetron{

  class EXPORTGADGETS_RADIAL gpuRadialSpiritPrepGadget : public gpuRadialPrepGadget
  {
    
  public:
    GADGET_DECLARE(gpuRadialSpiritPrepGadget);
    gpuRadialSpiritPrepGadget();
    virtual ~gpuRadialSpiritPrepGadget() {}
    
  protected:
    
    virtual int process_config(ACE_Message_Block *mb);

    virtual void reconfigure(unsigned int set, unsigned int slice, bool use_dcw = true );

    virtual boost::shared_ptr< hoNDArray<float_complext> > compute_csm( unsigned int buffer_idx );

    virtual boost::shared_ptr< hoNDArray<float_complext> > compute_reg( unsigned int set, 
                                                                        unsigned int slice, 
                                                                        bool new_frame );
    
    virtual void allocate_accumulation_buffer( unsigned int num_buffers );

    virtual cuBuffer<float,2>* get_buffer_ptr(int idx){
      return &this->acc_buffer_spirit_[idx];
    }
  };
}
