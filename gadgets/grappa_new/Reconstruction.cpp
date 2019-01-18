#include "Reconstruction.h"

#include "parallel/Merge.h"

namespace {
    using namespace Gadgetron::Core;
}

namespace Gadgetron::Grappa {
    GADGETRON_MERGE_EXPORT(Reconstruction);

    Reconstruction::Reconstruction(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : Merge(props) {}

    void Reconstruction::process(
            std::map<std::string, std::shared_ptr<Channel>> input,
            std::shared_ptr<Channel> output
    ) {
        TypedInputChannel<hoNDArray<std::complex<float>>> images{*input.at("images"), *output};
        TypedInputChannel<hoNDArray<float>> weights{*input.at("images"), *output};

        for (auto image : images) {
            GINFO_STREAM("Reconstruction handling image.")
        }
    }
}