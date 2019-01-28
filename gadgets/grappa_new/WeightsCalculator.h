#pragma once

#include "SliceAccumulator.h"

#include "Context.h"
#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Grappa {

    class WeightsCalculator : public Core::TypedGadgetNode<Slice> {
    public:
        WeightsCalculator(const Core::Context &, const std::unordered_map<std::string, std::string> &);

        NODE_PROPERTY(target_coils, int, "Number of target coils for GRAPPA recon", 0);

        // TODO: Get rid of these pointless things.
        NODE_PROPERTY(device_channels, int, "Number of device channels", 0);


        void process(Core::TypedInputChannel<Slice> &in, Core::OutputChannel &out) override;

    private:
        const Core::Context context;
    };
}
