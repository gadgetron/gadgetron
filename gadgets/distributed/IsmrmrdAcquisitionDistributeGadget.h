#ifndef ISMRMRDACQUISITIONDISTRIBUTEGADGET_H
#define ISMRMRDACQUISITIONDISTRIBUTEGADGET_H

#include "Gadget.h"
#include "gadgetron_distributed_gadgets_export.h"
#include "DistributeGadget.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class EXPORTDISTRIBUTEDGADGETS IsmrmrdAcquisitionDistributeGadget : public DistributeGadget
  {
  public:
    GADGET_DECLARE(IsmrmrdAcquisitionDistributeGadget);
    IsmrmrdAcquisitionDistributeGadget() {}
    virtual ~IsmrmrdAcquisitionDistributeGadget() {}

  protected:
    GADGET_PROPERTY_LIMITS(parallel_dimension, std::string,
      "Dimension that data will be parallelized over", "slice",
      GadgetPropertyLimitsEnumeration,
      "kspace_encode_step_1",
      "kspace_encode_step_2",
      "average",
      "slice",
      "contrast",
      "phase",
      "repetition",
      "set",
      "segment",
      "user_0",
      "user_1",
      "user_2",
      "user_3",
      "user_4",
      "user_5",
      "user_6",
      "user_7");


      virtual int node_index(ACE_Message_Block* m);
      virtual int message_id(ACE_Message_Block* m);

    };
  }
#endif //ISMRMRDACQUISITIONDISTRIBUTEGADGET_H
