#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <ismrmrd_hdf5.h>

namespace Gadgetron{

  class EXPORTGADGETSMRICORE CplxDumpGadget : 
  public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:

      CplxDumpGadget();
      ~CplxDumpGadget();

    protected:
      virtual int process_config(ACE_Message_Block* mb);

      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      virtual int close(unsigned long flags); //All the work is done here in this Gadget

    private:
      std::string filename_;
      ACE_Message_Queue<ACE_MT_SYNCH> buffer_;
    };
}
