#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <ctime>
#include "GadgetMRIHeaders.h"
#include "ismrmrd/meta.h"

namespace Gadgetron
{
    class NoiseSummaryGadget : public Core::ChannelGadget<void>
    {
    public:
        NoiseSummaryGadget(const Core::Context& context, const Core::GadgetProperties& props);
        
    	NODE_PROPERTY(noise_file, std::string, "Name of noise file", "");

        void process(Core::InputChannel<void>& input, Core::OutputChannel& output) override;

    private:
        Core::Context context;

    };
}
