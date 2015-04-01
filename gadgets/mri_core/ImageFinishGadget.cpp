#include "GadgetIsmrmrdReadWrite.h"
#include "ImageFinishGadget.h"

namespace Gadgetron{

    int ImageFinishGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
    {
        if (!this->controller_)
        {
            GERROR("Cannot return result to controller, no controller set");
            return -1;
        }

        GadgetContainerMessage<GadgetMessageIdentifier>* mb = new GadgetContainerMessage<GadgetMessageIdentifier>();

        mb->getObjectPtr()->id = GADGET_MESSAGE_ISMRMRD_IMAGE;
        mb->cont(m1);

        int ret = this->controller_->output_ready(mb);

        if ((ret < 0))
        {
            GERROR("Failed to return massage to controller\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(ImageFinishGadget);
}
