#pragma once

#include "SliceAccumulator.h"

#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Grappa {

    class ImageAccumulator : public Core::TypedGadgetNode<Slice> {
    public:
        ImageAccumulator(const Core::Context &, const std::unordered_map<std::string, std::string> &);

        void process(Core::TypedInputChannel<Slice> &in, Core::OutputChannel &out) override;

    private:
        const Core::Context context;
    };
}
