#include "ImageAccumulator.h"

#include "Node.h"

#include "log.h"

namespace {
    using namespace Gadgetron::Core;
}

namespace Gadgetron::Grappa {

    ImageAccumulator::ImageAccumulator(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode<Acquisition>(props) {}

    void ImageAccumulator::process(Core::TypedInputChannel<Core::Acquisition> &in, Core::OutputChannel &out) {
        GINFO_STREAM("Hello, I'm the ImageAccumulator process function. I'm running!");
    }

    GADGETRON_GADGET_EXPORT(ImageAccumulator);
}
