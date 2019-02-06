#pragma once

#include "SliceAccumulator.h"

#include "Context.h"
#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Grappa {

    class WeightsCalculator : public Core::TypedGadgetNode<Slice> {
    public:
        WeightsCalculator(const Core::Context &, const std::unordered_map<std::string, std::string> &);

        void process(Core::TypedInputChannel<Slice> &in, Core::OutputChannel &out) override;

    private:
        const Core::Context context;
    };
}
