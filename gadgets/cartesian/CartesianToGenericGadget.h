#ifndef CartesianToGenericGadget_H
#define CartesianToGenericGadget_H
#pragma once

#include "gadgetron_cartesian_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"

#include <ismrmrd.h>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>

namespace Gadgetron{

  class EXPORTGADGETS_CARTESIAN CartesianToGenericGadget :
    public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
  {

  public:
    GADGET_DECLARE(CartesianToGenericGadget);

    CartesianToGenericGadget();
    virtual ~CartesianToGenericGadget();

  protected:

    virtual int process_config(ACE_Message_Block* mb);
    
    virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);
    
  private:
    std::vector<unsigned int> matrix_size_;
    unsigned short center_phase_;
  };
}
#endif //CartesianToGenericGadget_H
