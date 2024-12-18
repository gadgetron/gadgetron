#pragma once

#include "Node.h"
#include "hoNDArray.h"

#include <complex>
#include <map>

namespace Gadgetron {



    class AcquisitionAccumulateTriggerGadget
        : public Core::ChannelGadget<std::variant<mrd::Acquisition, mrd::WaveformUint32>> {
    public:
        using Core::ChannelGadget<std::variant<mrd::Acquisition, mrd::WaveformUint32>>::ChannelGadget;
        void process(Core::InputChannel<std::variant<mrd::Acquisition, mrd::WaveformUint32>>& in,
            Core::OutputChannel& out) override;
        enum class TriggerDimension {
            kspace_encode_step_1,
            kspace_encode_step_2,
            average,
            slice,
            contrast,
            phase,
            repetition,
            set,
            segment,
            user_0,
            user_1,
            user_2,
            user_3,
            user_4,
            user_5,
            user_6,
            user_7,
            n_acquisitions,
            none
        };
        NODE_PROPERTY(trigger_dimension, TriggerDimension, "Dimension to trigger on", TriggerDimension::none);
        NODE_PROPERTY(sorting_dimension, TriggerDimension, "Dimension to trigger on", TriggerDimension::none);

        NODE_PROPERTY(n_acquisitions_before_trigger, unsigned long, "Number of acquisition before first trigger", 40);
        NODE_PROPERTY(n_acquisitions_before_ongoing_trigger, unsigned long, "Number of acquisition before ongoing triggers", 40);

        size_t trigger_events = 0;
    private:
        void send_data(Core::OutputChannel& out, std::map<unsigned int, mrd::AcquisitionBucket>& buckets,
                       std::vector<mrd::WaveformUint32>& waveforms);
    };

    void from_string(const std::string& str, AcquisitionAccumulateTriggerGadget::TriggerDimension& val);

}
