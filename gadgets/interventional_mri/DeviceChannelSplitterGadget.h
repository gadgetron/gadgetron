#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "gadgetron_interventional_mri_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  template <typename T> class EXPORTGADGETSINTERVENTIONAL_MRI DeviceChannelSplitterGadget : 
  public Gadget2<ISMRMRD::ImageHeader,hoNDArray< T > >
  {
  protected:
    virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, 
			GadgetContainerMessage< hoNDArray< T > >* m2);
  };
  
  class EXPORTGADGETSINTERVENTIONAL_MRI DeviceChannelSplitterGadgetUSHORT :
  public DeviceChannelSplitterGadget<uint16_t>
  {
  public:
    GADGET_DECLARE(DeviceChannelSplitterGadgetUSHORT);
  };

  class EXPORTGADGETSINTERVENTIONAL_MRI DeviceChannelSplitterGadgetFLOAT :
  public DeviceChannelSplitterGadget<float>
  {
  public:
    GADGET_DECLARE(DeviceChannelSplitterGadgetFLOAT);
  };

  class EXPORTGADGETSINTERVENTIONAL_MRI DeviceChannelSplitterGadgetCPLX :
  public DeviceChannelSplitterGadget< std::complex<float> >
  {
  public:
    GADGET_DECLARE(DeviceChannelSplitterGadgetCPLX);
  };
}

