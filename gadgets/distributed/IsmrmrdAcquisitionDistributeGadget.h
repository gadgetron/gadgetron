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
      
  protected:
    virtual int node_index(ACE_Message_Block* m);
    virtual int message_id(ACE_Message_Block* m);

  };
}
#endif //DISTRIBUTEGADGET_H
