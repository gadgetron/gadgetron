#ifndef CartesianToGenericGadget_H
#define CartesianToGenericGadget_H
#pragma once

#include "gadgetron_cartesian_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"

#include <ismrmrd/ismrmrd.h>
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
    GADGET_PROPERTY(matrix_size_as_multiple_of, int, "Force the matrix size to be a multiple of", 1);

    virtual int process_config(ACE_Message_Block* mb);    
    virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);
    
  private:
    std::vector<unsigned int> matrix_size_;
    unsigned short center_phase_;

    // We can enforce the encoding space dimension 
    // to be a multiple of the "warp size" (required for the gpu nfft)
    unsigned int warp_size_; 
  };
}
#endif //CartesianToGenericGadget_H
