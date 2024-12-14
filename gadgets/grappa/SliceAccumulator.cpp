#include "SliceAccumulator.h"

#include "common/AnnotatedAcquisition.h"
#include "common/grappa_common.h"

#include "Context.h"
#include "Channel.h"

#include "hoNDArray.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Grappa;
}

namespace Gadgetron::Grappa {

    SliceAccumulator::SliceAccumulator(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : ChannelGadget(context,props), context(context) {}

    void SliceAccumulator::process(Core::InputChannel<AnnotatedAcquisition> &in, Core::OutputChannel &out) {

        std::vector<AnnotatedAcquisition> acquisitions{};

        for (auto acquisition : in) {
            acquisitions.emplace_back(std::move(acquisition));

            if (is_last_in_slice(acquisitions.back())) {
                out.push(std::move(acquisitions));
                acquisitions = std::vector<AnnotatedAcquisition>{};
            }
        }
    }

    GADGETRON_GADGET_EXPORT(SliceAccumulator);
}