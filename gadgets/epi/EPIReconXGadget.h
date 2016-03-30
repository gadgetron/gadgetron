#ifndef EPIRECONXGADGET_H
#define EPIRECONXGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_epi_export.h"
#include "hoArmadillo.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

#include "EPIReconXObjectFlat.h"
#include "EPIReconXObjectTrapezoid.h"

namespace Gadgetron{

  class EXPORTGADGETS_EPI EPIReconXGadget : 
  public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      EPIReconXGadget();
      virtual ~EPIReconXGadget();
      
    protected:
      GADGET_PROPERTY(verboseMode, bool, "Verbose output", false);

      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      // in verbose mode, more info is printed out
      bool verboseMode_;

      // A set of reconstruction objects
      EPI::EPIReconXObjectTrapezoid<std::complex<float> > reconx;
      EPI::EPIReconXObjectFlat<std::complex<float> > reconx_other;

      // readout oversampling for reconx_other
      float oversamplng_ratio2_;

    };
}
#endif //EPIRECONXGADGET_H
