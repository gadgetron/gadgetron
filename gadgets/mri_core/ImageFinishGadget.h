/**
    \brief  Finalizes image and sends out generic message to generic channel
    \author Original: Thomas Sangild Sorensen
    \author PureGadget Conversion: David Christoffer Hansen
    \test   Tested by: simple_gre_3d.cfg and others
*/

#pragma once

#include "Gadget.h"
#include "Node.h"
#include "gadgetron_mricore_export.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE ImageFinishGadget : public Core::GenericChannelGadget {
    public:
        ImageFinishGadget(
                const Core::Context &context,
                const Core::GadgetProperties &properties
        ) : GenericChannelGadget(context,properties) {};

    protected:
        void process(Core::GenericInputChannel& in,
                    Core::OutputChannel& out) override;

    };
}

