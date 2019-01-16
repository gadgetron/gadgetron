#pragma once

#include <map>
#include <memory>

#include "parallel/Merge.h"
#include "Channel.h"

namespace Gadgetron::Grappa {
    class Reconstruction : public Core::Parallel::Merge {
    public:
        Reconstruction(const Core::Context &context, const std::unordered_map<std::string, std::string> &props);

        void process(
                std::map<std::string, std::shared_ptr<Core::Channel>> input,
                std::shared_ptr<Core::Channel> output
        ) override;
    };
}
