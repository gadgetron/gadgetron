#include "GadgetStreamController.h"
#include "gadgetheaders.h"
#include "GadgetStreamConfiguratorFactory.h"
#include "GadgetContainerMessage.h"
#include "GadgetSocketSender.h"
#include "NDArray.h"
#include "Gadget.h"

#include <complex>

int GadgetStreamController::open (void)
{
  ACE_TCHAR peer_name[MAXHOSTNAMELEN];
  ACE_INET_Addr peer_addr;
  if (this->sock_.get_remote_addr (peer_addr) == 0 &&
      peer_addr.addr_to_string (peer_name, MAXHOSTNAMELEN) == 0)
    ACE_DEBUG ((LM_DEBUG,
                ACE_TEXT ("(%P|%t) Connection from %s\n"),
                peer_name));

  if (!(output_ = new GadgetSocketSender(&this->sock_))) {
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("(%P|%t) %p\n"),
		       ACE_TEXT ("GadgetStreamController::open, Unable to create output_ task")), -1);

  }

  if (output_->register_writer(GADGET_MESSAGE_ACQUISITION, new GadgetAcquisitionMessageWriter()) < 0) {
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("(%P|%t) %p\n"),
		       ACE_TEXT ("GadgetStreamController::open, to register acquisition writer with output_ task")), -1);
  }

  if (output_->register_writer(GADGET_MESSAGE_IMAGE, new GadgetImageMessageWriter()) < 0) {
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("(%P|%t) %p\n"),
		       ACE_TEXT ("GadgetStreamController::open, to register acquisition writer with output_ task")), -1);
  }


  if (output_->open() == -1) {
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("(%P|%t) %p\n"),
		       ACE_TEXT ("GadgetStreamController::open, Unable to open output_ task")), -1);
  }

  return this->reactor ()->register_handler(this, 
					    ACE_Event_Handler::READ_MASK);  
}


int GadgetStreamController::handle_input (ACE_HANDLE)
{

  //Reading sequence:
  GadgetMessageIdentifier id;
  ssize_t recv_cnt = 0;
  if ((recv_cnt = this->sock_.recv_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
    ACE_DEBUG ((LM_DEBUG,
		ACE_TEXT ("(%P|%t) GadgetStreamController, unable to read message identifier\n")));
    return -1;
  }

  switch (id.id) {

  case GADGET_MESSAGE_ACQUISITION:
    if (read_acquisition() < 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("ACQ read failed\n")) );
      return -1;
    }
    break;
  case GADGET_MESSAGE_CONFIGURATION:
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Configuration received\n")) );
    if (read_configuration() < 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("CONFIG read failed\n")) );
      return -1;
    }
    break;
  case GADGET_MESSAGE_NEW_MEASUREMENT:
    break;
  case GADGET_MESSAGE_END_OF_SCAN:
    break;
  case GADGET_MESSAGE_IMAGE:
    break;
  case GADGET_MESSAGE_EMPTY:
    break;
  default:
    ACE_DEBUG( (LM_ERROR, ACE_TEXT("Received BAD MESSAGE IDENTIFIER %d\n"), id.id ) );
    return -1;
    break;
  }

  return 0;
}

int GadgetStreamController::handle_output (ACE_HANDLE)
{
  return 0;
}

int GadgetStreamController::handle_close (ACE_HANDLE, ACE_Reactor_Mask mask)
{
  if (mask == ACE_Event_Handler::WRITE_MASK)
    return 0;
  mask = ACE_Event_Handler::ALL_EVENTS_MASK |
         ACE_Event_Handler::DONT_CALL;
  this->reactor ()->remove_handler (this, mask);

  //We need to grab the pointer and set it to zero so that no output
  //is put on the output queue while shutting down.
  GadgetSocketSender* tmp = output_;
  output_ = 0;  
  tmp->close(1);
  delete tmp;

  this->stream_.close();
  this->sock_.close ();
  delete this;

  return 0;
}


int GadgetStreamController::read_configuration()
{
  GadgetMessageConfigurator c;
  ssize_t recv_cnt = 0;
  if ((recv_cnt = this->sock_.recv_n (&c, sizeof(GadgetMessageConfigurator))) <= 0) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) GadgetStreamController: Unable to read configuration\n")));
    return -1;
  }

  ACE_Message_Block* mb = new ACE_Message_Block(c.configuration_length);

  //ACE_TCHAR* config_info = 0;
  //ACE_NEW_RETURN(config_info, ACE_TCHAR[c.configuration_length],-1);

  if ((recv_cnt = this->sock_.recv_n (mb->wr_ptr(), c.configuration_length)) <= 0) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t)  GadgetStreamController: Unable to read configuration info\n")));
    return -1;
  }
  mb->wr_ptr(c.configuration_length);
  mb->set_flags(Gadget::GADGET_MESSAGE_CONFIG);

  //For now it is only possible to configure the stream once.
  if (!stream_configured_) {
    GadgetStreamConfigurator* cfg = 
      GadgetStreamConfiguratorFactory::CreateConfigurator(c,mb->rd_ptr(), this);
    
    auto_ptr<GadgetStreamConfigurator> co(cfg);
    if (cfg) {
      if (cfg->ConfigureStream(&this->stream_)) {
	mb->release();
	ACE_ERROR_RETURN( (LM_ERROR, ACE_TEXT("Unable to configure stream")), -1);
      }
    }
    co.release();
    
    stream_configured_ = true;
    if (cfg) delete cfg;
  }

  if (stream_.put(mb) < 0) {
    ACE_DEBUG( (LM_ERROR, 
		ACE_TEXT("Unable to send config down stream failed") ));
    mb->release();
    return -1;
  }

  return 0;
}

int GadgetStreamController::read_acquisition()
{
  
  GadgetContainerMessage<GadgetMessageAcquisition>* m1 = 
    new GadgetContainerMessage<GadgetMessageAcquisition>();

  GadgetContainerMessage<NDArray< std::complex<float> > >* m2 = 
    new GadgetContainerMessage< NDArray< std::complex<float> > >();

  m1->cont(m2);

  ssize_t recv_cnt = 0;
  if ((recv_cnt = 
       this->sock_.recv_n (m1->getObjectPtr(), 
			   sizeof(GadgetMessageAcquisition))) <= 0) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) Unable to read Acq header\n")));

    m1->release();

    return -1;
  }

  std::vector<int> adims;
  adims.push_back(m1->getObjectPtr()->samples);
  adims.push_back(m1->getObjectPtr()->channels);

  if (!m2->getObjectPtr()->create(adims)) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) Allocate sample data\n")));

    m1->release();

    return -1;
  }
  
  if ((recv_cnt = 
       this->sock_.recv_n
       (m2->getObjectPtr()->get_data_ptr(), 
	sizeof(std::complex<float>)*m1->getObjectPtr()->samples*m1->getObjectPtr()->channels)) <= 0) {

    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) Unable to read Acq data\n")));

    m1->release();

    return -1;
  }

  if (stream_.put(m1) < 0) {

    ACE_DEBUG( (LM_ERROR, 
		ACE_TEXT("Data send down stream failed") ));
    return -1;
  }

  return 0;
}

