#include "IsmrmrdAcquisitionDistributeGadget.h"
#include "GadgetMRIHeaders.h"

namespace Gadgetron{

  int IsmrmrdAcquisitionDistributeGadget::node_index(ACE_Message_Block* m)
  {
    auto h = AsContainerMessage<ISMRMRD::AcquisitionHeader>(m);

    if (!h) return GADGET_FAIL;

    std::string parallel_dimension_local = parallel_dimension.value();

    if (parallel_dimension_local.size() == 0) {
      return -1;
    } else if (parallel_dimension_local.compare("kspace_encode_step_1") == 0) {
      return h->getObjectPtr()->idx.kspace_encode_step_1;
    } else if (parallel_dimension_local.compare("kspace_encode_step_2") == 0) {
      return h->getObjectPtr()->idx.kspace_encode_step_2;
    } else if (parallel_dimension_local.compare("average") == 0) {
      return h->getObjectPtr()->idx.average;
    } else if (parallel_dimension_local.compare("slice") == 0) {
      return h->getObjectPtr()->idx.slice;
    } else if (parallel_dimension_local.compare("contrast") == 0) {
      return h->getObjectPtr()->idx.contrast;
    } else if (parallel_dimension_local.compare("phase") == 0) {
      return h->getObjectPtr()->idx.phase;
    } else if (parallel_dimension_local.compare("repetition") == 0) {
      return h->getObjectPtr()->idx.repetition;
    } else if (parallel_dimension_local.compare("set") == 0) {
      return h->getObjectPtr()->idx.set;
    } else if (parallel_dimension_local.compare("segment") == 0) {
      return h->getObjectPtr()->idx.segment;
    } else if (parallel_dimension_local.compare("user_0") == 0) {
      return h->getObjectPtr()->idx.user[0];
    } else if (parallel_dimension_local.compare("user_1") == 0) {
      return h->getObjectPtr()->idx.user[1];
    } else if (parallel_dimension_local.compare("user_2") == 0) {
      return h->getObjectPtr()->idx.user[2];
    } else if (parallel_dimension_local.compare("user_3") == 0) {
      return h->getObjectPtr()->idx.user[3];
    } else if (parallel_dimension_local.compare("user_4") == 0) {
      return h->getObjectPtr()->idx.user[4];
    } else if (parallel_dimension_local.compare("user_5") == 0) {
      return h->getObjectPtr()->idx.user[5];
    } else if (parallel_dimension_local.compare("user_6") == 0) {
      return h->getObjectPtr()->idx.user[6];
    } else if (parallel_dimension_local.compare("user_7") == 0) {
      return h->getObjectPtr()->idx.user[7];
    } else {
      GERROR("ERROR: Unkown parallel dimension\n");
      return -1;
    }
    return -1; //We should never reach this
  }

  int IsmrmrdAcquisitionDistributeGadget::message_id(ACE_Message_Block* m)
  {
    auto h = AsContainerMessage<ISMRMRD::AcquisitionHeader>(m);
    if (!h) return 0;

    return GADGET_MESSAGE_ISMRMRD_ACQUISITION;
  }

  GADGET_FACTORY_DECLARE(IsmrmrdAcquisitionDistributeGadget)

}
