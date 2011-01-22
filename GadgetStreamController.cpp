#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
#include "ace/OS_NS_stdio.h"
#include "ace/DLL.h"

#include <ticpp.h>

#include "GadgetStreamController.h"
#include "gadgetheaders.h"
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

  GadgetModule *head = 0;
  GadgetModule *tail = 0;

  if (tail == 0) {
    ACE_NEW_RETURN(tail, 
		   ACE_Module<ACE_MT_SYNCH>( ACE_TEXT("EndGadget"), 
					     new EndGadget() ),
		   -1);
    
    stream_.open(0,head,tail);
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
  case GADGET_MESSAGE_INITIALIZATION:
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Initialization received\n")) );
    if (read_initialization() < 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("INITIALIZATION read failed\n")) );
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

int GadgetStreamController::read_initialization()
{
  GadgetMessageInitializer ini;
  ssize_t recv_cnt = 0;
  if ((recv_cnt = this->sock_.recv_n (&ini, sizeof(GadgetMessageInitializer))) <= 0) {
    ACE_DEBUG ((LM_ERROR,
		ACE_TEXT ("(%P|%t) GadgetStreamController: Unable to read Initialization\n")));
    return GADGET_FAIL;
  }

  GADGET_DEBUG2("Received initialization file name: %s\n", ini.configuration_file);
  
  return this->configure(ini.configuration_file);
}

int GadgetStreamController::configure(char* init_filename)
{
  char * gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");

  if (!gadgetron_home) {
    GADGET_DEBUG1("GADGETRON_HOME environment variable not set\n");
    return GADGET_FAIL;
  }

  ACE_TCHAR config_file_name[4096];
  ACE_OS::sprintf(config_file_name, "%s/config/%s", gadgetron_home, init_filename);

  GADGET_DEBUG2("Running configuration: %s\n", config_file_name);

  try {
    ticpp::Document doc( config_file_name );
    doc.LoadFile();

    //Let's configure the stream
    ticpp::Iterator< ticpp::Node > child;
    ticpp::Element* pElem = doc.FirstChildElement("stream");
    for ( child = child.begin( pElem ); child != child.end(); child++ ) {
      if (ACE_OS::strncmp(child.Get()->ToElement()->Value().c_str(), "gadget", 6) == 0) {
	GADGET_DEBUG1("--Found gadget declaration\n");
	GADGET_DEBUG2("  Gadget Name: %s\n", child.Get()->ToElement()->GetAttribute("name").c_str());
	GADGET_DEBUG2("  Gadget dll: %s\n", child.Get()->ToElement()->GetAttribute("dll").c_str());
	GADGET_DEBUG2("  Gadget class: %s\n", child.Get()->ToElement()->GetAttribute("class").c_str());

	GadgetModule* m = create_gadget_module(child.Get()->ToElement()->GetAttribute("dll").c_str(),
					       child.Get()->ToElement()->GetAttribute("class").c_str(),
					       child.Get()->ToElement()->GetAttribute("name").c_str());

	if (!m) {
	  GADGET_DEBUG2("Failed to create GadgetModule from %s:%s\n", 
			child.Get()->ToElement()->GetAttribute("class").c_str(),
			child.Get()->ToElement()->GetAttribute("dll").c_str());
	  return GADGET_FAIL;
	}

	//Are there any attributes/parameters to set
	ticpp::Iterator< ticpp::Node > p;
	Gadget* g = dynamic_cast<Gadget*>(m->writer());
	for ( p = p.begin( child.Get()->ToElement() ); p != p.end(); p++ ) {
	  if (ACE_OS::strncmp(p.Get()->ToElement()->Value().c_str(), "property", 8) == 0) {
	    g->set_parameter(p.Get()->ToElement()->GetAttribute("name"),
			     p.Get()->ToElement()->GetAttribute("value"));
	  }
	}

	if (stream_.push(m) < 0) {
	  GADGET_DEBUG2("Failed to push Gadget %s onto stream\n", child.Get()->ToElement()->GetAttribute("name").c_str());
	}
      }
    }



  } catch( ticpp::Exception& ex ) {
    GADGET_DEBUG2("Error parsing XML config: %s\n", ex.what());
  }

  stream_configured_ = true;

  return GADGET_OK;
}

GadgetModule * GadgetStreamController::create_gadget_module(const char* DLL, 
							   const char* gadget, 
							   const char* gadget_module_name)
{

  ACE_DLL dll;

  ACE_TCHAR dllname[1024];
  ACE_OS::sprintf(dllname, "%s%s",ACE_DLL_PREFIX, DLL);

  ACE_TCHAR factoryname[1024];
  ACE_OS::sprintf(factoryname, "make_%s", gadget);

  int retval = dll.open (dllname);

  if (retval != 0)
    ACE_ERROR_RETURN ((LM_ERROR,
                       "%p, %s\n",
                       "dll.open", dllname),
                      0);

  //Function pointer
  typedef Gadget* (*GadgetCreator) (void);
  
  
  void *void_ptr = dll.symbol (factoryname);
  ptrdiff_t tmp = reinterpret_cast<ptrdiff_t> (void_ptr);
  GadgetCreator gc = reinterpret_cast<GadgetCreator> (tmp);

  if (gc == 0) {
    ACE_ERROR_RETURN ((LM_ERROR,
		       "%p, %s",
		       "dll.symbol", factoryname),
		      0);
  }

  dll.close ();

  Gadget* g = gc();
  
  if (!g) {
    GADGET_DEBUG1("Failed to create gadget using factory\n");
    return 0;
  }

  g->set_controller(this);

  GadgetModule *module = 0;
  ACE_NEW_RETURN (module,
		  GadgetModule (gadget_module_name, g), 
		  0);

  
  return module;
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
#if 0
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
#endif

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

