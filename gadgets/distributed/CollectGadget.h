#ifndef COLLECTGADGET_H
#define COLLECTGADGET_H

#include "Gadget.h"
#include "gadgetron_distributed_gadgets_export.h"

#include <complex>

namespace Gadgetron{

  class EXPORTDISTRIBUTEDGADGETS CollectGadget : public BasicPropertyGadget
  {
  public:
    GADGET_DECLARE(CollectGadget);

    CollectGadget();
    virtual ~CollectGadget();

  protected:
    GADGET_PROPERTY(pass_through_mode, bool,
      "If true, data will simply pass through to next gadget, otherwise return to controller", false);
    virtual int process(ACE_Message_Block* m);
    virtual int message_id(ACE_Message_Block* m);
  };
}
#endif //COLLECTGADGET_H
