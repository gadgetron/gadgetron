#ifndef ACCUMULATORGADGET_H
#define ACCUMULATORGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{
  
  class EXPORTGADGETSMRICORE AccumulatorGadget : 
  public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
      
    public:
      GADGET_DECLARE(AccumulatorGadget);
      
      AccumulatorGadget();
      ~AccumulatorGadget();
      
    protected:
      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);
      
      hoNDArray< std::complex<float> >* buffer_;
      std::vector<size_t> dimensions_;
      std::vector<float> field_of_view_;
      size_t slices_;
      long long image_counter_;
      long long image_series_;
    };
}
#endif //ACCUMULATORGADGET_H
