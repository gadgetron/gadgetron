#ifndef ACCUMULATORGADGET_H
#define ACCUMULATORGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{
  
  class EXPORTGADGETSMRICORE AccumulatorGadget : 
  public Gadget
    {
      
    public:
      GADGET_DECLARE(AccumulatorGadget);
      
      AccumulatorGadget();
      ~AccumulatorGadget();
      
    protected:

      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(ACE_Message_Block *mb);
      virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);
      
      boost::shared_ptr<hoNDArray< std::complex<float> > > buffer_;
      boost::shared_ptr< hoNDArray<float> > traj_buffer_;
      std::vector<size_t> trajectory_dimensions_;
      std::vector<size_t> dimensions_;
      std::vector<float> field_of_view_;
      size_t slices_;
      long long image_counter_;
      long long image_series_;
    };
}
#endif //ACCUMULATORGADGET_H
