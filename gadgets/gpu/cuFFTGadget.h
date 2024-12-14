#pragma once
#include "Gadget.h"
#include "hoNDArray.h"

#include <complex>

#include <boost/make_shared.hpp>


namespace Gadgetron{

  class cuFFTGadget :
    public Core::ChannelGadget<mrd::Image<std::complex<float>>>
    {
    public:
      using Core::ChannelGadget<mrd::Image<std::complex<float>>>::ChannelGadget;
      ~cuFFTGadget() override = default;
      void process(Core::InputChannel<mrd::Image<std::complex<float>>>& in, Core::OutputChannel& out) override;
  };
}
