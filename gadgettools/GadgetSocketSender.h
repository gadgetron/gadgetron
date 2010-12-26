#ifndef GADGETSOCKETSENDER_H
#define GADGETSOCKETSENDER_H

#include "ace/SOCK_Stream.h"
#include "ace/Task.h"

#include <complex>

#include "../gadgetheaders.h"
#include "NDArray.h"
#include "../GadgetContainerMessage.h"

/**
   Class for writing a specific message to a socket. 
   This is an abstract class, implementations need to be done for each message type.
 */
class GadgetMessageWriter
{
 public:
  /**
     Function must be implemented to read a specific message.
   */
  virtual int write(ACE_SOCK_Stream* stream, ACE_Message_Block* mb) = 0;
};

/**
   Default implementation of GadgetMessageWriter for Acquisition messages

 */
class GadgetAcquisitionMessageWriter : public GadgetMessageWriter
{

 public:
  virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb) 
  {
    GadgetContainerMessage<GadgetMessageAcquisition>* acqmb =
      dynamic_cast< GadgetContainerMessage<GadgetMessageAcquisition>* >(mb);
    
    GadgetContainerMessage< NDArray< std::complex<float> > >* datamb =
      dynamic_cast< GadgetContainerMessage< NDArray< std::complex<float> > >* >(acqmb->cont());
    
    if (!acqmb || !datamb) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), GadgetAcquisitionMessageWriter, invalid acquisition message objects")) );
       return -1;
    }

    ssize_t send_cnt = 0;
    
    GadgetMessageIdentifier id;
    id.id = GADGET_MESSAGE_ACQUISITION;
    
    if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
      ACE_DEBUG ((LM_ERROR,
		  ACE_TEXT ("(%P|%t) Unable to send acquisition message identifier\n")));
      
      return -1;
    }
    
    if ((send_cnt = sock->send_n (acqmb->getObjectPtr(), sizeof(GadgetMessageAcquisition))) <= 0) {
      ACE_DEBUG ((LM_ERROR,
		  ACE_TEXT ("(%P|%t) Unable to send acquisition header\n")));
      
      return -1;
    }
    
    if ((send_cnt = sock->send_n (datamb->getObjectPtr()->get_data_ptr(), sizeof(std::complex<float>)*datamb->getObjectPtr()->get_number_of_elements())) <= 0) {
      ACE_DEBUG ((LM_ERROR,
		  ACE_TEXT ("(%P|%t) Unable to send acquisition data\n")));
      
      return -1;
    }

    return 0;
  }

};

/**
   Default implementation of GadgetMessageWriter for Image messages

 */
class GadgetImageMessageWriter : public GadgetMessageWriter
{

 public:
  virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb) 
  {

    GadgetContainerMessage<GadgetMessageImage>* imagemb = 
      dynamic_cast< GadgetContainerMessage<GadgetMessageImage>* >(mb);
    
    GadgetContainerMessage< NDArray< std::complex<float> > >* datamb =
      dynamic_cast< GadgetContainerMessage< NDArray< std::complex<float> > >* >(imagemb->cont());
    
    if (!imagemb || !datamb) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), GadgetStreamController::handle_output, invalid image message objects")) );
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

};

class GadgetSocketSender : public ACE_Task<ACE_MT_SYNCH>
{

 public:
  typedef ACE_Task<ACE_MT_SYNCH> inherited;

  GadgetSocketSender(ACE_SOCK_Stream* socket) 
    : inherited()
    , socket_(socket)
  {
    ACE_TRACE(( ACE_TEXT("GadgetSocketSender::GadgetSocketSender") ));
  }

  virtual ~GadgetSocketSender() 
    {
      for (unsigned int i = 0; i < writers_.size(); i++) {
	if (writers_[i]) delete writers_[i]; 
      }  
    }

  virtual int init(void)
  {
    ACE_TRACE(( ACE_TEXT("Gadget::init") ));
    return 0;
  }

  virtual int open(void* = 0) 
  {
    ACE_TRACE(( ACE_TEXT("GagetSocketSender::open") ));
    
    return this->activate( THR_NEW_LWP | THR_JOINABLE,
			   1 );
  }


  virtual int register_writer(ACE_UINT16 slot, GadgetMessageWriter* writer) 
  {
    if (writer == 0) {
      return 0;
    }

    if (find_writer(slot) != 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GagetSocketSender, attempting to register writer in occupied slot\n")) );
      return -1;
    }
 
    slots_.push_back(slot);
    writers_.push_back(writer);
    return 0;
  }

  virtual int close(unsigned long flags)
  {
    ACE_TRACE(( ACE_TEXT("GagetSocketSender::close") ));
    return 0;
  }

  virtual int svc(void) 
  {
    ACE_TRACE(( ACE_TEXT("GadgetSocketSender::svc") ));

    ACE_Message_Block *mb;
    
    while (this->getq(mb) >= 0) {
      GadgetContainerMessage<GadgetMessageIdentifier>* mid;
      if (!(mid = dynamic_cast< GadgetContainerMessage<GadgetMessageIdentifier>* >(mb))) {
	ACE_ERROR_RETURN ((LM_ERROR,
			   ACE_TEXT ("(%P|%t) %p\n"),
			   ACE_TEXT ("GadgetSocketSender, Invalid message on output queue")), -1);
	
      }
    
      GadgetMessageWriter* w = find_writer(mid->getObjectPtr()->id);

      if (!w) {
	mb->release();
	ACE_ERROR_RETURN ((LM_ERROR,
			   ACE_TEXT ("(%P|%t) %p\n"),
			   ACE_TEXT ("GadgetSocketSender, Unable to find registered writer for message")), -1);
	
      }

      if (w->write(socket_, mb->cont()) < 0) {
	mb->release();
	ACE_ERROR_RETURN ((LM_ERROR,
			   ACE_TEXT ("(%P|%t) %p\n"),
			   ACE_TEXT ("GadgetSocketSender, failed to write message to socket")), -1);
      } 

      mb->release();
    }
    return 0;
  }


 protected:
  ACE_SOCK_Stream* socket_;
  std::vector<ACE_UINT16> slots_;
  std::vector<GadgetMessageWriter*> writers_;

  GadgetMessageWriter* find_writer(ACE_UINT16 slot) {
    GadgetMessageWriter* ret = 0;
    for (unsigned int i = 0; i < slots_.size(); i++) {
      if (slots_[i] == slot) ret = writers_[i]; break;
    }
    return ret;
  }
};

#endif //GADGETSOCKETSENDER_H
