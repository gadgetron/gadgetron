/**
    \brief  Passes through an acquisition to the next gadget in the pipeline
    \test   Tested by: simple_gre_acquisition_passthrough.cfg
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"
#include "gadgetron_debugging_export.h"


namespace Gadgetron{
  class AcquisitionPassthroughGadget : public Core::ChannelGadget<Core::Acquisition> 
    {
      public:
        using Core::ChannelGadget<Core::Acquisition>::ChannelGadget;
        ~AcquisitionPassthroughGadget() override = default;
        void process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) override;
    };
}