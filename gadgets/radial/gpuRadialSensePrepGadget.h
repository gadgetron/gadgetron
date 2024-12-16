#pragma once

#include "gpuRadialPrepGadget.h"

namespace Gadgetron{

  class gpuRadialSensePrepGadget : public gpuRadialPrepGadget
  {

  public:
    gpuRadialSensePrepGadget() : gpuRadialPrepGadget() {}
    virtual ~gpuRadialSensePrepGadget() {}

  protected:

    virtual void reconfigure(unsigned int set, unsigned int slice, bool use_dcw = true);

    virtual boost::shared_ptr< hoNDArray<float_complext> > compute_csm( unsigned int buffer_idx );

    virtual boost::shared_ptr< hoNDArray<float_complext> > compute_reg( unsigned int set,
                                                                        unsigned int slice,
                                                                        bool new_frame );

    virtual void allocate_accumulation_buffer( unsigned int num_buffers );

    virtual cuBuffer<float,2>* get_buffer_ptr(int idx){
      return (this->buffer_using_solver_) ? &this->acc_buffer_sense_cg_[idx] : &this->acc_buffer_sense_[idx];
    }
  };
}
