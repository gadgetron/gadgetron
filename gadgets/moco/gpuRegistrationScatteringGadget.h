#ifndef gpuRegistrationScatteringGadget_H
#define gpuRegistrationScatteringGadget_H

#include "cuNDArray_operators.h"
#include "cuNDArray_utils.h"
#include "cuCKOpticalFlowSolver.h"
#include "RegistrationScatteringGadget.h"

namespace Gadgetron{  

  class EXPORTGADGETS_MOCO gpuRegistrationScatteringGadget2D :
    public RegistrationScatteringGadget< cuNDArray<float>, 2 >
  {    
  public:
    GADGET_DECLARE(gpuRegistrationScatteringGadget2D);
    gpuRegistrationScatteringGadget2D() : RegistrationScatteringGadget< cuNDArray<float>, 2 >() {}
    virtual ~gpuRegistrationScatteringGadget2D() {}

  protected:
    virtual int setup_solver();
    virtual int set_continuation( GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, cuNDArray<float> *continuation );
    virtual int write_displacement_field( cuNDArray<float> *displacements );
  };
}

#endif //gpuRegistrationScatteringGadget_H
