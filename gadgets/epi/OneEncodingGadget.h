/** \file   OneEncodingGadget.h
    \brief  This is the class gadget to make sure EPI Flash Ref lines are in the same encoding space as the imaging lines.
    \author Hui Xue
*/

#ifndef ONEENCODINGGADGET_H
#define ONEENCODINGGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"

#include <complex>

namespace Gadgetron {

    class OneEncodingGadget :
        public Gadget1<mrd::Acquisition>
    {
    public:
        OneEncodingGadget();
        virtual ~OneEncodingGadget();

    protected:
        virtual int process(GadgetContainerMessage< mrd::Acquisition>* m1);
    };
}

#endif // ONEENCODINGGADGET_H
