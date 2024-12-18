#include "AcquisitionWaveformBranch.h"

using namespace Gadgetron::Core::Parallel;

namespace {
    std::string select_channel(const mrd::Acquisition &) { return "acquisitions"; }
    std::string select_channel(const mrd::WaveformUint32 &) { return "waveforms"; }
}

namespace Gadgetron::Examples {

    AcquisitionWaveformBranch::AcquisitionWaveformBranch(
            const Core::Context &,
            const Core::GadgetProperties &properties
    ) : TypedBranch<AcquisitionOrWaveform>(properties) {}

    void AcquisitionWaveformBranch::process(
        Core::InputChannel<AcquisitionOrWaveform> &input,
            std::map<std::string, Core::OutputChannel> output
    ) {
        for (auto acq_or_wav : input) {
            auto &channel = output.at(visit([](auto &aw) { return select_channel(aw); }, acq_or_wav));
            channel.push(std::move(acq_or_wav));
        }
    }

    GADGETRON_BRANCH_EXPORT(AcquisitionWaveformBranch)
}