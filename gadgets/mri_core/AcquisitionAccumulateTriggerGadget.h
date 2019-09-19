#pragma once

#include "Node.h"
#include "hoNDArray.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <map>
#include "mri_core_acquisition_bucket.h"

namespace Gadgetron{

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
    void from_string(const std::string& str, TriggerDimension& val);



  class AcquisitionAccumulateTriggerGadget :
public Core::ChannelGadget<Core::variant<Core::Acquisition ,Core::Waveform>>
    {
    public:
        using Core::ChannelGadget<Core::variant<Core::Acquisition, Core::Waveform>>::ChannelGadget;
        void process(Core::InputChannel<Core::variant<Core::Acquisition, Core::Waveform>>& in,
            Core::OutputChannel& out) override;

      NODE_PROPERTY(trigger_dimension, TriggerDimension, "Dimension to trigger on", TriggerDimension::none);
    NODE_PROPERTY(sorting_dimension, TriggerDimension, "Dimension to trigger on", TriggerDimension::none);

      NODE_PROPERTY(n_acquisitions_before_trigger, unsigned long, "Number of acquisition before first trigger", 40);
      NODE_PROPERTY(n_acquisitions_before_ongoing_trigger, unsigned long, "Number of acquisition before ongoing triggers", 40);
      


    };

  
}
