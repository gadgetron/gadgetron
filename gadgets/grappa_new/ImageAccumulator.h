#pragma once

#include <chrono>

#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Grappa {

    class ImageAccumulator : public Core::TypedGadgetNode<Core::Acquisition> {
    public:
        ImageAccumulator(const Core::Context &, const std::unordered_map<std::string, std::string> &);

        void process(Core::TypedInputChannel<Core::Acquisition> &in, Core::OutputChannel &out) override;

    private:
        const Core::Context context;
    };
}
