#include "GadgetStreamController.h"
#include "gadgetheaders.h"
#include "GadgetStreamConfiguratorFactory.h"
#include "GadgetContainerMessage.h"
#include "NDArray.h"

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

  return this->reactor ()->register_handler(this, ACE_Event_Handler::READ_MASK);
}


int GadgetStreamController::handle_input (ACE_HANDLE)
{

  //Reading sequence:
  GadgetMessageIdentifier id;
  ssize_t recv_cnt = 0;
  if ((recv_cnt = this->sock_.recv_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
    ACE_DEBUG ((LM_DEBUG,
		ACE_TEXT ("(%P|%t) GadgetStreamController, unable to read message identifier\n")));
    stream_.close();
    return -1;
  }

  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("GadgetMessageIdentifier received: %d\n"), id.id) );

  switch (id.id) {

  case GADGET_MESSAGE_ACQUISITION:
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Acquisition received\n")) );
    if (read_acquisition() < 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("ACQ read failed\n")) );
      stream_.close();
      return -1;
    }
    break;
  case GADGET_MESSAGE_CONFIGURATION:
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Configuration received\n")) );
    if (read_configuration() < 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("CONFIG read failed\n")) );
      stream_.close();
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
    ACE_DEBUG( (LM_ERROR, ACE_TEXT("Received BAD MESSAGE IDENTIFIER\n")) );
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
  this->sock_.close ();
  delete this;
  return 0;
}


int GadgetStreamController::read_configuration()
{
  GadgetMessageConfigurator c;
  ssize_t recv_cnt = 0;
  if ((recv_cnt = this->sock_.recv_n (&c, sizeof(GadgetMessageConfigurator))) <= 0) {
    ACE_DEBUG ((LM_DEBUG,
		ACE_TEXT ("(%P|%t) GadgetStreamController: Unable to read configuration\n")));
    return -1;
  }

  ACE_TCHAR* config_info = 0;
  ACE_NEW_RETURN(config_info, ACE_TCHAR[c.configuration_length],-1);

  if ((recv_cnt = this->sock_.recv_n (config_info, c.configuration_length)) <= 0) {
    ACE_DEBUG ((LM_DEBUG,
		ACE_TEXT ("(%P|%t)  GadgetStreamController: Unable to read configuration info\n")));
    return -1;
  }


  GadgetStreamConfigurator* cfg = 
    GadgetStreamConfiguratorFactory::CreateConfigurator(c,config_info);

  auto_ptr<GadgetStreamConfigurator> co(cfg);
  if (cfg) {
    if (cfg->ConfigureStream(&this->stream_)) {
      delete [] config_info;
      ACE_ERROR_RETURN( (LM_ERROR, ACE_TEXT("Unable to configure stream")), -1);
    }
  }
  co.release();
  
  if (cfg) delete cfg;
  delete [] config_info;

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
    ACE_DEBUG ((LM_DEBUG,
		ACE_TEXT ("(%P|%t) Unable to read Acq header\n")));

    m1->release();

    return -1;
  }

  std::cout << "Samples: " << m1->getObjectPtr()->samples << std::endl;

  std::vector<int> adims;
  adims.push_back(m1->getObjectPtr()->samples);

  if (!m2->getObjectPtr()->create(adims)) {
    ACE_DEBUG ((LM_DEBUG,
		ACE_TEXT ("(%P|%t) Allocate sample data\n")));

    m1->release();

    return -1;
  }
  
  if ((recv_cnt = 
       this->sock_.recv_n
       (m2->getObjectPtr()->get_data_ptr(), 
	sizeof(std::complex<float>)*m1->getObjectPtr()->samples)) <= 0) {

    ACE_DEBUG ((LM_DEBUG,
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
