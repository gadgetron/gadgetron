#include "SliceAccumulator.h"

#include "common/slice.h"

#include "Context.h"
#include "Channel.h"
#include "Types.h"

#include "hoNDArray.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;
}

namespace Gadgetron::Grappa {

    SliceAccumulator::SliceAccumulator(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode(props), context(context) {}

    void SliceAccumulator::process(TypedInputChannel<Acquisition> &in, OutputChannel &out) {

        std::vector<Acquisition> acquisitions{};

        for (auto acquisition : in) {
            acquisitions.emplace_back(std::move(acquisition));

            if (is_last_in_slice(acquisitions.back())) {
                GINFO_STREAM("Finished accumulating slice: " << slice_of(acquisitions.back()));

                out.push(std::move(acquisitions));
                acquisitions = std::vector<Acquisition>{};

            }
        }
    }

    GADGETRON_GADGET_EXPORT(SliceAccumulator);
}