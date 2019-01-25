#include "Combine.h"

#include "parallel/Merge.h"

namespace {
    using namespace Gadgetron::Core;
}

namespace Gadgetron::Grappa {
    GADGETRON_MERGE_EXPORT(Combine);

    Combine::Combine(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : Merge(props) {}

    void Combine::process(
            std::map<std::string, std::shared_ptr<Channel>> input,
            std::shared_ptr<Channel> output
    ) {
        TypedInputChannel<Image>   images{*input.at("images"), *output};
        TypedInputChannel<Weights> weights{*input.at("images"), *output};

        for (auto image : images) {
            GINFO_STREAM("Reconstruction handling image.")

            // Excellent!
        }
    }
}