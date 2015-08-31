#include "IsmrmrdImageDistributeGadget.h"
#include "GadgetMRIHeaders.h"

namespace Gadgetron{

  int IsmrmrdImageDistributeGadget::node_index(ACE_Message_Block* m)
  {
    auto h = AsContainerMessage<ISMRMRD::ImageHeader>(m);

    if (!h) return GADGET_FAIL;

    std::string parallel_dimension_local = parallel_dimension.value();

    if (parallel_dimension_local.size() == 0) {
      return -1;
    } else if (parallel_dimension_local.compare("average") == 0) {
      return h->getObjectPtr()->average;
    } else if (parallel_dimension_local.compare("slice") == 0) {
      return h->getObjectPtr()->slice;
    } else if (parallel_dimension_local.compare("contrast") == 0) {
      return h->getObjectPtr()->contrast;
    } else if (parallel_dimension_local.compare("phase") == 0) {
      return h->getObjectPtr()->phase;
    } else if (parallel_dimension_local.compare("repetition") == 0) {
      return h->getObjectPtr()->repetition;
    } else if (parallel_dimension_local.compare("set") == 0) {
      return h->getObjectPtr()->set;
    } else {
      GERROR("ERROR: Unkown parallel dimension\n");
      return -1;
    }
    return -1; //We should never reach this
  }

  int IsmrmrdImageDistributeGadget::message_id(ACE_Message_Block* m)
  {
    auto h = AsContainerMessage<ISMRMRD::ImageHeader>(m);
    if (!h) return 0;

    return GADGET_MESSAGE_ISMRMRD_IMAGE;
  }

  GADGET_FACTORY_DECLARE(IsmrmrdImageDistributeGadget)

}
