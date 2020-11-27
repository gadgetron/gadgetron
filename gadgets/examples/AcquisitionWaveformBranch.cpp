#include "AcquisitionWaveformBranch.h"

using namespace Gadgetron::Core;
using namespace Gadgetron::Core::Parallel;

namespace {
    std::string select_channel(const Acquisition &) { return "acquisitions"; }
    std::string select_channel(const Waveform &) { return "waveforms"; }
}

namespace Gadgetron::Examples {

    AcquisitionWaveformBranch::AcquisitionWaveformBranch(
            const Context &,
            const GadgetProperties &properties
    ) : TypedBranch<AcquisitionOrWaveform>(properties) {}

    void AcquisitionWaveformBranch::process(
        InputChannel<AcquisitionOrWaveform> &input,
            std::map<std::string, OutputChannel> output
    ) {
        for (auto acq_or_wav : input) {
            auto &channel = output.at(visit([](auto &aw) { return select_channel(aw); }, acq_or_wav));
            channel.push(std::move(acq_or_wav));
        }
    }

    GADGETRON_BRANCH_EXPORT(AcquisitionWaveformBranch)
}