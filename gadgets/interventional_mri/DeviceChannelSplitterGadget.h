#pragma once

#include "Gadget.h"
#include "hoNDArray.h"

#include <complex>

namespace Gadgetron{

  template <typename T> class DeviceChannelSplitterGadget :
  public Gadget1<mrd::Image<T>>
  {
  protected:
    virtual int process(GadgetContainerMessage<mrd::Image<T>>* m1);
  };

  class DeviceChannelSplitterGadgetUSHORT :
  public DeviceChannelSplitterGadget<uint16_t>
  {
  };

  class DeviceChannelSplitterGadgetFLOAT :
  public DeviceChannelSplitterGadget<float>
  {
  };

  class DeviceChannelSplitterGadgetCPLX :
  public DeviceChannelSplitterGadget< std::complex<float> >
  {
  };
}

