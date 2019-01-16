#include "WeightsCalculator.h"

#include "Gadget.h"

namespace {
    using namespace Gadgetron::Core;
}

namespace Gadgetron::Grappa {

    WeightsCalculator::WeightsCalculator(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode<Acquisition>(props) {}

    void WeightsCalculator::process(Core::TypedInputChannel<Core::Acquisition> &in, Core::OutputChannel &out) {
        GINFO_STREAM("Hello, I'm the WeightsCalculator process function. I'm running!");
    }

    GADGETRON_GADGET_EXPORT(WeightsCalculator);
}