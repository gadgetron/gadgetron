#include "GadgetIsmrmrdReadWrite.h"
#include "ImageFinishAttribGadget.h"

namespace Gadgetron
{
    template <typename T>
    int ImageFinishAttribGadget<T>::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< T > >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3)
    {
        if (!this->controller_)
        {
            GERROR("Cannot return result to controller, no controller set");
            return -1;
        }

        GadgetContainerMessage<GadgetMessageIdentifier>* mb = new GadgetContainerMessage<GadgetMessageIdentifier>();

        switch (sizeof(T))
        {
        case 2: //Unsigned short
        {
            if (typeid(T) == typeid(unsigned short))
                mb->getObjectPtr()->id = GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT;
            else
                mb->getObjectPtr()->id = GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_SHORT;
        }
            break;
        case 4: //Float
            mb->getObjectPtr()->id = GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT;
            break;
        case 8: //Complex float
            mb->getObjectPtr()->id = GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT;
            break;
        default:
            GDEBUG("Wrong data size detected: %d\n", sizeof(T));
            mb->release();
            m1->release();
            return GADGET_FAIL;
        }

        mb->cont(m1);

        int ret = this->controller_->output_ready(mb);

        if ((ret < 0))
        {
            GDEBUG("Failed to return massage to controller\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    //Declare factories for the various template instances
    GADGET_FACTORY_DECLARE(ImageFinishAttribGadgetFLOAT)
    GADGET_FACTORY_DECLARE(ImageFinishAttribGadgetUSHORT)
    GADGET_FACTORY_DECLARE(ImageFinishAttribGadgetSHORT)
    GADGET_FACTORY_DECLARE(ImageFinishAttribGadgetCPLX)
}
