#include "AcquisitionPassthroughGadget.h"

namespace Gadgetron{
  void Gadgetron::AcquisitionPassthroughGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) {
      for (auto acquisition : in) {
          out.push(acquisition);
      }  
  }
  GADGETRON_GADGET_EXPORT(AcquisitionPassthroughGadget)
}