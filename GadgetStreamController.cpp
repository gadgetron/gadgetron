#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
#include "ace/OS_NS_stdio.h"
#include "ace/DLL.h"

#include <ticpp.h>

#include "GadgetStreamController.h"
#include "GadgetContainerMessage.h"
#include "NDArray.h"
#include "Gadget.h"
#include "GadgetronRuntimeLinking.h"

#include <complex>

int GadgetStreamController::open (void)
{

  //We will set up the controllers message queue such that when a packet is enqueued write will be triggered.
  this->notifier_.reactor (this->reactor ());
  this->msg_queue ()->notification_strategy (&this->notifier_);

  ACE_TCHAR peer_name[MAXHOSTNAMELEN];
  ACE_INET_Addr peer_addr;
  if (peer().get_remote_addr (peer_addr) == 0 &&
      peer_addr.addr_to_string (peer_name, MAXHOSTNAMELEN) == 0)
    ACE_DEBUG ((LM_DEBUG,
                ACE_TEXT ("(%P|%t) Connection from %s\n"),
                peer_name));

  //We have to have these basic types to be able to receive configuration file for stream
  readers_.insert(GADGET_MESSAGE_CONFIG_FILE, 
		  new GadgetMessageConfigFileReader());

  readers_.insert(GADGET_MESSAGE_PARAMETER_SCRIPT, 
		  new GadgetMessageScriptReader());


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
					    ACE_Event_Handler::READ_MASK);// | ACE_Event_Handler::WRITE_MASK);  
}


int GadgetStreamController::handle_input (ACE_HANDLE)
{

  //Reading sequence:
  GadgetMessageIdentifier id;
  ssize_t recv_cnt = 0;
  if ((recv_cnt = peer().recv_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
    ACE_DEBUG ((LM_DEBUG,
		ACE_TEXT ("(%P|%t) GadgetStreamController, unable to read message identifier\n")));
    return -1;
  }

  GadgetMessageReader* r = readers_.find(id.id);

  if (!r) {
    GADGET_DEBUG2("Unrecognized Message ID received: %d\n", id.id);
    return GADGET_FAIL;
  }

  ACE_Message_Block* mb = r->read(&peer());

  if (!mb) {
    GADGET_DEBUG1("GadgetMessageReader returned null pointer\n");
    return GADGET_FAIL;
  }

  //We need to handle some special cases to make sure that we can get a stream set up.
  if (id.id == GADGET_MESSAGE_CONFIG_FILE) {
    GadgetContainerMessage<GadgetMessageConfigurationFile>* cfgm = 
      AsContainerMessage<GadgetMessageConfigurationFile>(mb);

    if (!cfgm) {
      GADGET_DEBUG1("Failed to cast message block to configuration file\n");
      mb->release();
      return GADGET_FAIL;
    } else {
      if (this->configure(cfgm->getObjectPtr()->configuration_file) != GADGET_OK) {
	GADGET_DEBUG1("GadgetStream configuration failed\n");
	mb->release();
	return GADGET_FAIL;
      } else {
	return GADGET_OK;
      }
    }

  }
    
 
  ACE_Time_Value wait = ACE_OS::gettimeofday() + ACE_Time_Value(0,10000); //10ms from now
  if (stream_.put(mb) == -1) {
    GADGET_DEBUG2("Failed to put stuff on stream, too long wait, %d\n",  ACE_OS::last_error () ==  EWOULDBLOCK);
    mb->release();
    return GADGET_FAIL;
  } 

  return GADGET_OK;

}


int GadgetStreamController::output_ready(ACE_Message_Block* mb) { 
  
  int res = this->putq(mb);
  return res;
}

int GadgetStreamController::handle_output (ACE_HANDLE)
{
  ACE_Message_Block *mb = 0;
  ACE_Time_Value nowait (ACE_OS::gettimeofday ());

  //Send packages as long as we have them
  while (-1 != this->getq (mb, &nowait)) {
    GadgetContainerMessage<GadgetMessageIdentifier>* mid =
      AsContainerMessage<GadgetMessageIdentifier>(mb);

    if (!mid) {
      GADGET_DEBUG1("Invalid message on output queue\n");
      mb->release();
      return GADGET_FAIL;
    }

    GadgetMessageWriter* w = writers_.find(mid->getObjectPtr()->id);
    
    if (!w) {
      GADGET_DEBUG2("Unrecognized Message ID received: %d\n", mid->getObjectPtr()->id);
      return GADGET_FAIL;
    }
    
    if (w->write(&peer(),mb->cont()) < 0) {
      GADGET_DEBUG1("Failed to write Message using writer\n");
      mb->release ();
      return GADGET_FAIL;
    }

    mb->release();
  }


  if (this->msg_queue ()->is_empty ()) {
    //No point in coming back to handle_ouput until something is put on the queue,
    //in which case, the msg queue's notification strategy will tell us
    this->reactor ()->cancel_wakeup(this, ACE_Event_Handler::WRITE_MASK);
  } else {
    //There is still more on the queue, let's come back when idle
    this->reactor ()->schedule_wakeup(this, ACE_Event_Handler::WRITE_MASK);
  }

  return 0;
}

int GadgetStreamController::handle_close (ACE_HANDLE, ACE_Reactor_Mask mask)
{

  GADGET_DEBUG1("handle_close called\n");

  if (mask == ACE_Event_Handler::WRITE_MASK)
    return 0;

  GADGET_DEBUG1("Shutting down stream and closing up shop\n");

  mask = ACE_Event_Handler::ALL_EVENTS_MASK |
         ACE_Event_Handler::DONT_CALL;
  this->reactor ()->remove_handler (this, mask);


  this->stream_.close();
  delete this;
  return 0;
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

    //Configuration of readers
    ticpp::Iterator< ticpp::Node > child;
    ticpp::Element* pElem = doc.FirstChildElement("readers", false);
    if (pElem) {
      for ( child = child.begin( pElem ); child != child.end(); child++ ) {
	if (child.Get()->Type() == TiXmlNode::ELEMENT &&
	    ACE_OS::strncmp(child.Get()->ToElement()->Value().c_str(), "reader", 6) == 0) {
	  GADGET_DEBUG1("--Found reader declaration\n");
	  GADGET_DEBUG2("  Reader dll: %s\n", child.Get()->ToElement()->GetAttribute("dll").c_str());
	  GADGET_DEBUG2("  Reader class: %s\n", child.Get()->ToElement()->GetAttribute("class").c_str());
	  GADGET_DEBUG2("  Reader slot: %d\n", ACE_OS::atoi(child.Get()->ToElement()->GetAttribute("slot").c_str()));

	  GadgetMessageReader* r =
	    load_dll_component<GadgetMessageReader>(child.Get()->ToElement()->GetAttribute("dll").c_str(),
							child.Get()->ToElement()->GetAttribute("class").c_str());
       
	  if (!r) {
	    GADGET_DEBUG1("Failed to load GadgetMessageReader from DLL\n");
	    return GADGET_FAIL;
	  }

	  readers_.insert(ACE_OS::atoi(child.Get()->ToElement()->GetAttribute("slot").c_str()),
			  r);

	}
      }
    }
    //Configuration of readers end
    

    //Configuration of writers
    //Configuration of readers
    pElem = doc.FirstChildElement("writers", false);
    if (pElem) {
      for ( child = child.begin( pElem ); child != child.end(); child++ ) {
	if (child.Get()->Type() == TiXmlNode::ELEMENT &&
	    ACE_OS::strncmp(child.Get()->ToElement()->Value().c_str(), "writer", 6) == 0) {
	  GADGET_DEBUG1("--Found writer declaration\n");
	  GADGET_DEBUG2("  Writer dll: %s\n", child.Get()->ToElement()->GetAttribute("dll").c_str());
	  GADGET_DEBUG2("  Writer class: %s\n", child.Get()->ToElement()->GetAttribute("class").c_str());
	  GADGET_DEBUG2("  Writer slot: %d\n", ACE_OS::atoi(child.Get()->ToElement()->GetAttribute("slot").c_str()));

	  GadgetMessageWriter* w =
	    load_dll_component<GadgetMessageWriter>(child.Get()->ToElement()->GetAttribute("dll").c_str(),
							child.Get()->ToElement()->GetAttribute("class").c_str());
       
	  if (!w) {
	    GADGET_DEBUG1("Failed to load GadgetMessageWriter from DLL\n");
	    return GADGET_FAIL;
	  }

	  writers_.insert(ACE_OS::atoi(child.Get()->ToElement()->GetAttribute("slot").c_str()),
			  w);

	}
      }
    }

    //Configuration of writers end

    //Let's configure the stream
    pElem = doc.FirstChildElement("stream");
    for ( child = child.begin( pElem ); child != child.end(); child++ ) {
      if (child.Get()->Type() == TiXmlNode::ELEMENT && ACE_OS::strncmp(child.Get()->ToElement()->Value().c_str(), "gadget", 6) == 0) {
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

  Gadget* g = load_dll_component<Gadget>(DLL,gadget);
  
  if (!g) {
    GADGET_DEBUG1("Failed to load gadget using factory\n");
    return 0;
  }

  g->set_controller(this);

  GadgetModule *module = 0;
  ACE_NEW_RETURN (module,
		  GadgetModule (gadget_module_name, g), 
		  0);

  
  return module;
}


template <class T>  
T* GadgetStreamController::load_dll_component(const char* DLL, const char* component_name)
{

  ACE_DLL_Manager* dllmgr = ACE_DLL_Manager::instance();

  ACE_DLL_Handle* dll;
  ACE_SHLIB_HANDLE dll_handle = 0;

  ACE_TCHAR dllname[1024];
  ACE_OS::sprintf(dllname, "%s%s",ACE_DLL_PREFIX, DLL);

  ACE_TCHAR factoryname[1024];
  ACE_OS::sprintf(factoryname, "make_%s", component_name);


  dll = dllmgr->open_dll (dllname, ACE_DEFAULT_SHLIB_MODE, dll_handle );

  if (!dll)
    ACE_ERROR_RETURN ((LM_ERROR,
                       "%p, ---%s---\n",
                       "dll.open", dllname),
                      0);

  //Function pointer
  typedef T* (*ComponentCreator) (void);
  
  
  void *void_ptr = dll->symbol (factoryname);
  ptrdiff_t tmp = reinterpret_cast<ptrdiff_t> (void_ptr);
  ComponentCreator cc = reinterpret_cast<ComponentCreator> (tmp);

  if (cc == 0) {
    ACE_ERROR_RETURN ((LM_ERROR,
		       "%p,  ---%s---\n",
		       "dll.symbol", factoryname),
		      0);
  }


  T* c = cc();
  
  if (!c) {
    GADGET_DEBUG1("Failed to create component using factory\n");
    return 0;
  }

  return c;
}
