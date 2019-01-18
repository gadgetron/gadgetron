#include "ImageAccumulator.h"

#include <chrono>

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

    void ImageAccumulator::process(TypedInputChannel<Acquisition> &in, OutputChannel &out) {

        GINFO_STREAM("Hello, I'm the ImageAccumulator process function. I'm running!");

        std::vector<hoNDArray<std::complex<float>>> image_data;
        std::vector<std::chrono::milliseconds> time_stamps;

        for (auto acquisition : in) {
            GINFO_STREAM("ImageAccumulator handling an acquisition.")
        }

        GINFO_STREAM("ImageAccumulator acquisition loop over - channel closed.")
    }

    GADGETRON_GADGET_EXPORT(ImageAccumulator);
}
