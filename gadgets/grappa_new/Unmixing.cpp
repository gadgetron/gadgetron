#include "Unmixing.h"

#include "parallel/Merge.h"

namespace {
    using namespace Gadgetron::Core;
}

namespace Gadgetron::Grappa {
    GADGETRON_MERGE_EXPORT(Unmixing);

    Unmixing::Unmixing(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : Merge(props) {}

    void Unmixing::process(
            std::map<std::string, InputChannel> input,
            OutputChannel output
    ) {
        TypedInputChannel<Image> images{input.at("images"), output};
        TypedInputChannel<Weights> weights{input.at("weights"), output};

        for (auto image : images) {
            GINFO_STREAM("Received unmixing job. Excellent.");
        }
    }
}