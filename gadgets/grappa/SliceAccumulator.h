#pragma once

#include <vector>

#include "ChannelReorderer.h"

#include "Node.h"

namespace Gadgetron::Grappa {

    using Slice = std::vector<AnnotatedAcquisition>;

    class SliceAccumulator : public Core::ChannelGadget<AnnotatedAcquisition> {
    public:
        SliceAccumulator(const Core::Context &, const std::unordered_map<std::string, std::string> &);

        void process(Core::InputChannel<AnnotatedAcquisition> &in, Core::OutputChannel &out) override;

    private:
        const Core::Context context;
    };
}

