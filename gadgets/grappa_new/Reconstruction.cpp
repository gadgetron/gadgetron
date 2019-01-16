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

    }
}