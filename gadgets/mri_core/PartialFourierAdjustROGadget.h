
#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"

namespace Gadgetron
{

/// for incoming readout
/// if not the noise scan and the partial fourier along readout is detected
/// the readout data will be realigned with center of echo at the centre of incoming 1D array
class PartialFourierAdjustROGadget : public Gadgetron::Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
{
public:

    PartialFourierAdjustROGadget();

protected:

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(Gadgetron::GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
        Gadgetron::GadgetContainerMessage< Gadgetron::hoNDArray< std::complex<float> > >* m2);

    unsigned int maxRO_;
};

}
