#pragma once

#include <vector>

#include "Types.h"
#include "Node.h"

namespace Gadgetron::Grappa {

    using Slice = std::vector<Core::Acquisition>;

    class SliceAccumulator : public Core::TypedGadgetNode<Core::Acquisition> {
    public:
        SliceAccumulator(const Core::Context &, const std::unordered_map<std::string, std::string> &);

        void process(Core::TypedInputChannel<Core::Acquisition> &in, Core::OutputChannel &out) override;

    private:
        const Core::Context context;
    };
}

