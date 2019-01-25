#include "WeightsCalculator.h"

#include "common/AcquisitionBuffer.h"
#include "Combine.h"

#include "Gadget.h"

namespace {
    using namespace Gadgetron::Core;
}

namespace Gadgetron::Grappa {

    WeightsCalculator::WeightsCalculator(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode<Acquisition>(props), context(context) {}

    void WeightsCalculator::process(Core::TypedInputChannel<Core::Acquisition> &in, Core::OutputChannel &out) {

        AcquisitionBuffer buffer{context};

        auto tuple = in.try_pop();

        // All right guys! Here we go! We need to somehow just sort shit out. We need to do some things:
        // Read until the input queue is empty - try_pop should.

        // Read acquisitions until channel is empty.
        // Keep an eye on the acuisitions - determine if we discard weights buffer (?), or recalculate.

        // recalculate if needed; emit new weights.
        // Wait for more acquisitions.



        for (const auto &acquisition : in) {
            buffer.add_acquisition(acquisition);
        }
    }

    GADGETRON_GADGET_EXPORT(WeightsCalculator);
}