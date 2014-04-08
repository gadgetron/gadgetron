#ifndef EPIRECONXGADGET_H
#define EPIRECONXGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"
#include "hoArmadillo.h"

#include <ismrmrd.h>
#include <complex>

#include "EPIReconXObjectTrapezoid.h"

namespace Gadgetron{

  class EXPORTGADGETSMRICORE EPIReconXGadget : 
  public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      EPIReconXGadget();
      virtual ~EPIReconXGadget();
      
    protected:
      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      // in verbose mode, more info is printed out
      bool verboseMode_;

      // A set of reconstruction objects
      EPI::EPIReconXObjectTrapezoid<std::complex<float> > reconx;
      //std::vector< EPI::EPIReconXObjectTrapezoid<std::complex<float> > > reconx;
    };
}
#endif //EPIRECONXGADGET_H
