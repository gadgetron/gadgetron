/**
    \brief  
    \author Original: Thomas Sangild Sorensen
    \author PureGadget Conversion: Kristoffer Langeland Knudsen
    \test   Untested
*/

#pragma once
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace Gadgetron{
  
    class FlowPhaseSubtractionGadget : public Core::ChannelGadget<Core::Image<std::complex<float>>>
    {

    public:
        using Core::ChannelGadget<Core::Image<std::complex<float>>>::ChannelGadget;

        ~FlowPhaseSubtractionGadget() override = default;

        void process(Core::InputChannel<Core::Image<std::complex<float>>>& in, Core::OutputChannel& out) override;
    };
}

