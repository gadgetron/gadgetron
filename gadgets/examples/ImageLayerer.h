#pragma once

#include "parallel/Merge.h"
#include "Gadget.h"
#include "Types.h"

namespace Gadgetron::Examples {

    class ImageLayerer : public Core::Parallel::Merge {
    public:
        ImageLayerer(const Core::Context &, const Core::GadgetProperties &);
        void process(std::map<std::string, Core::GenericInputChannel>, Core::OutputChannel) override;
    };
}
