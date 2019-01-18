#pragma once

#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Grappa {

    class WeightsCalculator : public Core::TypedGadgetNode<Core::Acquisition> {
    public:
        WeightsCalculator(const Core::Context &, const std::unordered_map<std::string, std::string> &);

        void process(Core::TypedInputChannel<Core::Acquisition> &in, Core::OutputChannel &out) override;
    };
}
