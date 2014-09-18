#ifndef gpuRegistrationAveragingGadget_H
#define gpuRegistrationAveragingGadget_H

#include "cuNDArray_operators.h"
#include "cuNDArray_utils.h"
#include "cuCKOpticalFlowSolver.h"
#include "RegistrationAveragingGadget.h"

namespace Gadgetron{  

  class EXPORTGADGETS_MOCO gpuRegistrationAveragingGadget2D :
    public RegistrationAveragingGadget< cuNDArray<float>, 2 >
  {    
  public:
    GADGET_DECLARE(gpuRegistrationAveragingGadget2D);

    gpuRegistrationAveragingGadget2D() : RegistrationAveragingGadget< cuNDArray<float>, 2 >() {}
    virtual ~gpuRegistrationAveragingGadget2D() {}

  protected:
    virtual int setup_solver();
    virtual int set_continuation( GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, cuNDArray<float> *continuation );
  };
}

#endif //gpuRegistrationAveragingGadget_H
