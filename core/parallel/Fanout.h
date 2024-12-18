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
        void process(InputChannel<ARGS...> &, std::map<std::string, OutputChannel>) override;
    };

    using AcquisitionFanout = Core::Parallel::Fanout<mrd::Acquisition>;
    using WaveformFanout = Core::Parallel::Fanout<mrd::WaveformUint32>;
    using ImageFanout = Core::Parallel::Fanout<mrd::AnyImage>;
}

#include "Fanout.hpp"
