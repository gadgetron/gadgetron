#pragma once

#include "parallel/Branch.h"

namespace Gadgetron::Examples {

    using AcquisitionOrWaveform = Core::variant<Core::Acquisition, Core::Waveform>;

    class AcquisitionWaveformBranch : public Core::Parallel::TypedBranch<AcquisitionOrWaveform> {
    public:
        AcquisitionWaveformBranch(const Core::Context &, const Core::GadgetProperties &);
        void process(
                Core::TypedInputChannel<AcquisitionOrWaveform> &,
                std::map<std::string, Core::OutputChannel>
        ) override;
    };
}