#include "SliceAccumulator.h"

#include "Context.h"
#include "Channel.h"
#include "Types.h"

#include "hoNDArray.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;

    bool is_last_in_slice(const Acquisition &acquisition) {
        auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        return header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
    }
}

namespace Gadgetron::Grappa {

    SliceAccumulator::SliceAccumulator(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode(props), context(context) {}

    void SliceAccumulator::process(TypedInputChannel<Acquisition> &in, OutputChannel &out) {

        std::vector<Acquisition> acquisitions;

        for (auto acquisition : in) {
            auto last_in_slice = is_last_in_slice(acquisition);
            acquisitions.emplace_back(std::move(acquisition));

            if (last_in_slice) {
                out.push(std::move(acquisitions));
                acquisitions = std::vector<Acquisition>{};
            }
        }
    }

    GADGETRON_GADGET_EXPORT(SliceAccumulator);
}