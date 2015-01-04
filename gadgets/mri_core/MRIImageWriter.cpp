#include "GadgetIsmrmrdReadWrite.h"
#include "MRIImageWriter.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "log.h"

#include <complex>

namespace Gadgetron{

  template <typename T>
  int MRIImageWriter<T>::write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb)
  {
    GadgetContainerMessage<ISMRMRD::ImageHeader>* imagemb =
      AsContainerMessage<ISMRMRD::ImageHeader>(mb);

    if (!imagemb) {
      GERROR("MRIImageWriter::write, invalid image message objects, 1\n");
      return -1;
    }

    GadgetContainerMessage< hoNDArray< T > >* datamb =
      AsContainerMessage< hoNDArray< T > >(imagemb->cont());

    if (!datamb) {
      GERROR("MRIImageWriter::write, invalid image message objects\n");
      return -1;
    }

    ssize_t send_cnt = 0;
    GadgetMessageIdentifier id;
    switch (sizeof(T)) {
    case 2: //Unsigned short
      id.id = GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT;
      break;
    case 4: //Float
      id.id = GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT;
      break;
    case 8: //Complex float
      id.id = GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT;
      break;
    default:
      GERROR("MRIImageWriter Wrong data size detected\n");
      return GADGET_FAIL;
    }


    //Let's check if the image header is consistent with the data array size before sending:
    unsigned long expected_elements = imagemb->getObjectPtr()->matrix_size[0]*
      imagemb->getObjectPtr()->matrix_size[1] *  imagemb->getObjectPtr()->matrix_size[2]*imagemb->getObjectPtr()->channels;

    if (expected_elements !=  datamb->getObjectPtr()->get_number_of_elements()) {
      GDEBUG("Number of header elements %d is inconsistent with number of elements in NDArray %d\n",expected_elements, datamb->getObjectPtr()->get_number_of_elements());
      GDEBUG("Header dimensions: %d, %d, %d, %d\n",imagemb->getObjectPtr()->matrix_size[0],imagemb->getObjectPtr()->matrix_size[1],imagemb->getObjectPtr()->matrix_size[2], imagemb->getObjectPtr()->channels);
      GDEBUG("Number of array dimensions: %d:\n", datamb->getObjectPtr()->get_number_of_dimensions());
      for (size_t i = 0; i < datamb->getObjectPtr()->get_number_of_dimensions(); i++) {
	GDEBUG("Dimensions %d: %d\n", i, datamb->getObjectPtr()->get_size(i));
      }
      return -1;
    }

    if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
      GERROR("Unable to send image message identifier\n");
      return -1;
    }

    if ((send_cnt = sock->send_n ( imagemb->getObjectPtr(), sizeof(ISMRMRD::ImageHeader))) <= 0) {
      GERROR("Unable to send image header\n");
      return -1;
    }

    if ((send_cnt = sock->send_n (datamb->getObjectPtr()->get_data_ptr(), sizeof(T)*datamb->getObjectPtr()->get_number_of_elements())) <= 0) {
      GERROR("Unable to send image data\n");
      return -1;
    }

    return 0;
  }

  GADGETRON_WRITER_FACTORY_DECLARE(MRIImageWriterFLOAT);
  GADGETRON_WRITER_FACTORY_DECLARE(MRIImageWriterUSHORT);
  GADGETRON_WRITER_FACTORY_DECLARE(MRIImageWriterCPLX);
}
