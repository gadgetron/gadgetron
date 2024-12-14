#pragma once

#include "parallel/Branch.h"

namespace Gadgetron::Examples {

    using AcquisitionOrWaveform = std::variant<mrd::Acquisition, mrd::WaveformUint32>;

    class AcquisitionWaveformBranch : public Core::Parallel::TypedBranch<AcquisitionOrWaveform> {
    public:
        AcquisitionWaveformBranch(const Core::Context &, const Core::GadgetProperties &);
        void process(
                Core::InputChannel<AcquisitionOrWaveform> &,
                std::map<std::string, Core::OutputChannel>
        ) override;
    };
}