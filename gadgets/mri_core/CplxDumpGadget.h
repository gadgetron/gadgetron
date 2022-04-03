/** 
    \brief  Dumps complex acquisition data to file
    \author Original: Thomas Sangild
    \author PureGadget Conversion: Kristoffer Langeland Knudsen
    \test   Untested 
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"
#include <ismrmrd/ismrmrd.h>

namespace Gadgetron{

  class CplxDumpGadget :
    public Core::ChannelGadget<Core::Acquisition>
    {
    public:
      using Core::ChannelGadget<Core::Acquisition>::ChannelGadget;
      ~CplxDumpGadget() override = default;

      NODE_PROPERTY(filename, std::string, "Filename of dumpfile", "profiles.cplx");
      void process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) override;
    };
}
