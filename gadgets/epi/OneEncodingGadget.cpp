#include "OneEncodingGadget.h"

namespace Gadgetron {

    OneEncodingGadget::OneEncodingGadget()
    {
    }

    OneEncodingGadget::~OneEncodingGadget()
    {
    }

    int OneEncodingGadget::process(GadgetContainerMessage< mrd::Acquisition>* m1)
    {
        m1->getObjectPtr()->head.encoding_space_ref = 0;

        if (this->next()->putq(m1) < 0)
        {
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(OneEncodingGadget)
}
