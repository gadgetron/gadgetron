
#include "Gadget.h"
#include "GadgetMessageInterface.h"
#include "GadgetStreamInterface.h"

#include "EndGadget.h"

namespace Gadgetron {

    int EndGadget::close(unsigned long flags)
    {
        GDEBUG("Close called in EndGadget with flags %d\n", flags);

        GadgetContainerMessage<GadgetMessageIdentifier>* mb =
                new GadgetContainerMessage<GadgetMessageIdentifier>();

        mb->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

        if (controller_->output_ready(mb) < 0) {
            return GADGET_FAIL;
        }

        GDEBUG("Calling close in base class  with flags %d\n", flags);
        return Gadget::close(flags);
    }

    int EndGadget::process(ACE_Message_Block *m)
    {
        m->release();
        return 0;
    }

    int EndGadget::next_step(ACE_Message_Block *m)
    {
        m->release();
        return 0;
    }

    int EndGadget::process_config(ACE_Message_Block * m) {
        m->release();
        return 0;
    }
}
