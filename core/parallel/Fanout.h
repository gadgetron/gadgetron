#pragma once

#include <map>
#include <memory>

#include "Branch.h"

#include "Channel.h"

namespace Gadgetron::Core::Parallel {

    template<class ...ARGS>
    class Fanout : public TypedBranch<ARGS...> {
    public:
        Fanout(const Context &context, const GadgetProperties &props);
        void process(TypedInputChannel<ARGS...> &, std::map<std::string, OutputChannel>) override;
    };
}

#include "Fanout.hpp"
