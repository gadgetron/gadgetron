#include "OneEncodingGadget.h"

namespace Gadgetron {

    OneEncodingGadget::OneEncodingGadget()
    {
    }

    OneEncodingGadget::~OneEncodingGadget()
    {
    }

    int OneEncodingGadget::process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
    {
        m1->getObjectPtr()->encoding_space_ref = 0;

        if (this->next()->putq(m1) < 0)
        {
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(OneEncodingGadget)
}
