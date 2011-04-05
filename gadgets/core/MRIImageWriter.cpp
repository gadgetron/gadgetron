#include "MRIImageWriter.h"

#include <complex>

#include "GadgetContainerMessage.h"
#include "hoNDArray.h"


int MRIImageWriter::write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb) 
{
  
  GadgetContainerMessage<GadgetMessageImage>* imagemb = 
    AsContainerMessage<GadgetMessageImage>(mb);
  
  if (!imagemb) {
    ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), MRIImageWriter::write, invalid image message objects, 1\n")) );
    return -1;    
  }

  GadgetContainerMessage< hoNDArray< std::complex<float> > >* datamb =
    AsContainerMessage< hoNDArray< std::complex<float> > >(imagemb->cont());
  
  if (!datamb) {
    ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), MRIImageWriter::write, invalid image message objects\n")) );
    return -1;
  }
  
  
  ssize_t send_cnt = 0;
  GadgetMessageIdentifier id;
  id.id = GADGET_MESSAGE_IMAGE;
  
  if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) Unable to send image message identifier\n")));
    
    return -1;
  }
  
  if ((send_cnt = sock->send_n (imagemb->getObjectPtr(), sizeof(GadgetMessageImage))) <= 0) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) Unable to send image header\n")));
      
    return -1;
  }

  if ((send_cnt = sock->send_n (datamb->getObjectPtr()->get_data_ptr(), sizeof(std::complex<float>)*datamb->getObjectPtr()->get_number_of_elements())) <= 0) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) Unable to send image data\n")));
      
    return -1;
  }

  return 0;
}

GADGETRON_WRITER_FACTORY_DECLARE(MRIImageWriter)
