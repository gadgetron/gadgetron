#ifndef EPIRECONXGADGET_H
#define EPIRECONXGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "hoArmadillo.h"

#include <complex>

#include "EPIReconXObjectFlat.h"
#include "EPIReconXObjectTrapezoid.h"

namespace Gadgetron{

  class EPIReconXGadget : 
  public Gadget1<mrd::Acquisition>
    {
    public:
      EPIReconXGadget();
      virtual ~EPIReconXGadget();
      
    protected:
      GADGET_PROPERTY(verboseMode, bool, "Verbose output", false);

      virtual int process_config(const mrd::Header& header);
      virtual int process(GadgetContainerMessage<mrd::Acquisition>* m1);

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
