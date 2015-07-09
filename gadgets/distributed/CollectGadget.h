#ifndef COLLECTGADGET_H
#define COLLECTGADGET_H

#include "Gadget.h"
#include "gadgetron_distributed_gadgets_export.h"

#include <complex>

namespace Gadgetron{

  class EXPORTDISTRIBUTEDGADGETS CollectGadget : public Gadget
    {
    public:
      GADGET_DECLARE(DistributeGadget);
      
    protected:
      virtual int process(ACE_Message_Block* m);
    };
}
#endif //DISTRIBUTEGADGET_H
