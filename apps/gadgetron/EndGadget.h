#ifndef ENDGADGET_H_
#define ENDGADGET_H_

#include "Gadget.h"
#include "GadgetMessageInterface.h"

namespace Gadgetron{
    class EXPORTGADGETBASE EndGadget : public Gadget
    {
        virtual int close(unsigned long flags);

    protected:
        virtual int process(ACE_Message_Block *m);
        virtual int next_step(ACE_Message_Block *m);
        virtual int process_config(ACE_Message_Block * m);
    };
}

#endif /* ENDGADGET_H_ */
