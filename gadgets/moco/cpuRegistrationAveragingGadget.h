#ifndef cpuRegistrationAveragingGadget_H
#define cpuRegistrationAveragingGadget_H

#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"
#include "hoCKOpticalFlowSolver.h"
#include "RegistrationAveragingGadget.h"

namespace Gadgetron{  

  class EXPORTGADGETS_MOCO cpuRegistrationAveragingGadget2D :
    public RegistrationAveragingGadget< hoNDArray<float>, 2 >
  {    
  public:
    GADGET_DECLARE(cpuRegistrationAveragingGadget2D);
    cpuRegistrationAveragingGadget2D() : RegistrationAveragingGadget< hoNDArray<float>, 2 >() {}
    virtual ~cpuRegistrationAveragingGadget2D() {}

  protected:
    virtual int setup_solver();
    virtual int set_continuation( GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, hoNDArray<float> *continuation );
  };
}

#endif //cpuRegistrationAveragingGadget_H
