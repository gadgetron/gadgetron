#ifndef ISMRMRDIMAGEDISTRIBUTEGADGET_H
#define ISMRMRDIMAGEDISTRIBUTEGADGET_H

#include "Gadget.h"
#include "gadgetron_distributed_gadgets_export.h"
#include "DistributeGadget.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class EXPORTDISTRIBUTEDGADGETS IsmrmrdImageDistributeGadget : public DistributeGadget
  {
  public:
    GADGET_DECLARE(IsmrmrdImageDistributeGadget);
    IsmrmrdImageDistributeGadget() {}
    virtual ~IsmrmrdImageDistributeGadget() {}

  protected:
    GADGET_PROPERTY_LIMITS(parallel_dimension, std::string,
      "Dimension that data will be parallelized over", "slice",
      GadgetPropertyLimitsEnumeration,
      "average",
      "slice",
      "contrast",
      "phase",
      "repetition",
      "set");

    virtual int node_index(ACE_Message_Block* m);
    virtual int message_id(ACE_Message_Block* m);
  };
}
#endif //ISMRMRDIMAGEDISTRIBUTEGADGET_H
