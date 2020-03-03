#pragma once

#include <Context.h>
#include <PropertyMixin.h>
#include <Channel.h>

#include "Merge.h"

namespace Gadgetron::Core::Parallel {

    class UnorderedMerge : public Merge {
    public:
        UnorderedMerge(const Context &context, const GadgetProperties &props);
        void process(std::map<std::string, GenericInputChannel>, OutputChannel) override;
    };
}