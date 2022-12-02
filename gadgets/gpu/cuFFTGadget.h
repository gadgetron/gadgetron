#pragma once
#include "Gadget.h"
#include "hoNDArray.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

#include <boost/make_shared.hpp>


namespace Gadgetron{

  class cuFFTGadget :
    public Core::ChannelGadget<Core::Image<std::complex<float>>>
    {
    public:
      using Core::ChannelGadget<Core::Image<std::complex<float>>>::ChannelGadget;
      ~cuFFTGadget() override = default;
      void process(Core::InputChannel<Core::Image<std::complex<float>>>& in, Core::OutputChannel& out) override;
  };
}
