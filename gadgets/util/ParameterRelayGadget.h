#ifndef PARAMETERRELAYGADGET_H
#define PARAMETERRELAYGADGET_H

#include "Gadget.h"
#include "gadgetron_util_gadgets_export.h"

namespace Gadgetron{

  class EXPORTUTILGADGETS ParameterRelayGadget : public Gadget
    {
    public:
      GADGET_DECLARE(ParameterRelayGadget);
      
    protected:
      virtual int process(ACE_Message_Block* m);
    };
}
#endif //PARAMETERRELAYGADGET_H
