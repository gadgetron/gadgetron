#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
#include "ace/OS_NS_stdio.h"
#include "ace/DLL.h"
#include "ace/DLL_Manager.h"
#include "ace/OS_NS_netdb.h"

#include <tinyxml.h>

#include "GadgetStreamController.h"
#include "GadgetContainerMessage.h"
#include "Gadget.h"

#include <complex>
#include <fstream>

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

  readers_.insert(GADGET_MESSAGE_CONFIG_SCRIPT, 
		  new GadgetMessageScriptReader());

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
      if (this->configure_from_file(std::string(cfgm->getObjectPtr()->configuration_file)) != GADGET_OK) {
	GADGET_DEBUG1("GadgetStream configuration failed\n");
	mb->release();
	return GADGET_FAIL;
      } else {
	mb->release();
	return GADGET_OK;
      }
    }
  } else if (id.id == GADGET_MESSAGE_CONFIG_SCRIPT) {
    std::string xml_config(mb->rd_ptr(), mb->length());
    if (this->configure(xml_config) != GADGET_OK) {
      GADGET_DEBUG1("GadgetStream configuration failed\n");
      mb->release();
      return GADGET_FAIL;
    } else {
      mb->release();
      return GADGET_OK;
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

int GadgetStreamController::configure_from_file(std::string config_xml_filename)
{

  char * gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");
  ACE_TCHAR config_file_name[4096];
  ACE_OS::sprintf(config_file_name, "%s/config/%s", gadgetron_home, config_xml_filename.c_str());

  GADGET_DEBUG2("Running configuration: %s\n", config_file_name);
  
  std::ifstream file (config_file_name, std::ios::in|std::ios::binary|std::ios::ate);
  if (file.is_open())
  {
    size_t size = file.tellg();
    char* buffer = new char [size];
    if (!buffer) {
      GADGET_DEBUG1("Unable to create temporary buffer for configuration file\n");
      return GADGET_FAIL;
    }
    file.seekg (0, std::ios::beg);
    file.read (buffer, size);
    file.close();
    std::string xml_file_contents(buffer,size);
    delete[] buffer;

    return configure(xml_file_contents);

  } else {
    GADGET_DEBUG2("Unable to open configuation file: %s\n", config_file_name);
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

int GadgetStreamController::configure(std::string config_xml_string)
{
  /*
  char * gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");

  if (!gadgetron_home) {
    GADGET_DEBUG1("GADGETRON_HOME environment variable not set\n");
    return GADGET_FAIL;
  }

  ACE_TCHAR config_file_name[4096];
  ACE_OS::sprintf(config_file_name, "%s/config/%s", gadgetron_home, init_filename);

  GADGET_DEBUG2("Running configuration: %s\n", config_file_name);
  */

  TiXmlDocument doc;//( config_file_name );
  doc.Parse(config_xml_string.c_str());

  //Configuration of readers
  TiXmlNode* child = 0;
  TiXmlNode* pElem = doc.FirstChild("readers");
  if (pElem) {
    while( (child = pElem->IterateChildren(child)) ) {
      if ((child->Type() == TiXmlNode::TINYXML_ELEMENT) &&
	  (ACE_OS::strncmp(child->ToElement()->Value(),"reader",6) == 0)) 
	{
	  std::string dllname(child->ToElement()->Attribute("dll"));
	  std::string classname(child->ToElement()->Attribute("class"));
	  int slot (ACE_OS::atoi(child->ToElement()->Attribute("slot")));

	  GADGET_DEBUG1("--Found reader declaration\n");
	  GADGET_DEBUG2("  Reader dll: %s\n", dllname.c_str());
	  GADGET_DEBUG2("  Reader class: %s\n", classname.c_str());
	  GADGET_DEBUG2("  Reader slot: %d\n", slot);

	  GadgetMessageReader* r =
	    load_dll_component<GadgetMessageReader>(dllname.c_str(),
						    classname.c_str());
       
	  if (!r) {
	    GADGET_DEBUG1("Failed to load GadgetMessageReader from DLL\n");
	    return GADGET_FAIL;
	  }

	  readers_.insert(slot, r);
	}
    }
  } 
  //Configuration of readers end
    

  //Configuration of writers
  pElem = doc.FirstChild("writers");
  child = 0;
  if (pElem) {
    while( (child = pElem->IterateChildren(child)) ) {
      if ((child->Type() == TiXmlNode::TINYXML_ELEMENT) &&
	  (ACE_OS::strncmp(child->ToElement()->Value(),"writer",6) == 0)) 
	{
	  std::string dllname(child->ToElement()->Attribute("dll"));
	  std::string classname(child->ToElement()->Attribute("class"));
	  int slot (ACE_OS::atoi(child->ToElement()->Attribute("slot")));
	  
	  GADGET_DEBUG1("--Found writer declaration\n");
	  GADGET_DEBUG2("  Writer dll: %s\n", dllname.c_str());
	  GADGET_DEBUG2("  Writer class: %s\n", classname.c_str());
	  GADGET_DEBUG2("  Writer slot: %d\n", slot);

	  GadgetMessageWriter* w =
	    load_dll_component<GadgetMessageWriter>(dllname.c_str(),
						    classname.c_str());
       
	  if (!w) {
	    GADGET_DEBUG1("Failed to load GadgetMessageWriter from DLL\n");
	    return GADGET_FAIL;
	  }

	  writers_.insert(slot, w);
	}
    }
  } 
  //Configuration of writers end

  //Let's configure the stream
  pElem = doc.FirstChild("stream");
  child = 0;
  if (pElem) {
    while( (child = pElem->IterateChildren(child)) ) {
      if ((child->Type() == TiXmlNode::TINYXML_ELEMENT) &&
	  (ACE_OS::strncmp(child->ToElement()->Value(),"gadget",6) == 0)) 
	{
	  std::string gadgetname(child->ToElement()->Attribute("name"));
	  std::string dllname(child->ToElement()->Attribute("dll"));
	  std::string classname(child->ToElement()->Attribute("class"));
	  
	  GADGET_DEBUG1("--Found gadget declaration\n");
	  GADGET_DEBUG2("  Gadget Name: %s\n", gadgetname.c_str());
	  GADGET_DEBUG2("  Gadget dll: %s\n", dllname.c_str());
	  GADGET_DEBUG2("  Gadget class: %s\n", classname.c_str());
	  
	  GadgetModule* m = create_gadget_module(dllname.c_str(),
						 classname.c_str(),
						 gadgetname.c_str());
	  
	  if (!m) {
	    GADGET_DEBUG2("Failed to create GadgetModule from %s:%s\n", 
			  classname.c_str(),
			  dllname.c_str());
	    return GADGET_FAIL;
	  }

	  TiXmlNode* child_param = 0;

	  Gadget* g = dynamic_cast<Gadget*>(m->writer());//Get the gadget out of the module

	  while( (child_param = child->IterateChildren(child_param)) ) {
	    if (child_param->Type() == TiXmlNode::TINYXML_ELEMENT &&
		(ACE_OS::strncmp(child_param->ToElement()->Value(), "property", 8) == 0)) 
	      {
		std::string propertyname(child_param->ToElement()->Attribute("name"));
		std::string propertyvalue(child_param->ToElement()->Attribute("value"));
		
		g->set_parameter(propertyname.c_str(), propertyvalue.c_str());
	      }
	    
	  }

	  if (stream_.push(m) < 0) {
	    GADGET_DEBUG2("Failed to push Gadget %s onto stream\n", gadgetname.c_str());
	    delete m;
	    return GADGET_FAIL;
	  }
	  
	}
    }
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

  ACE_DLL_Handle* dll = 0;
  ACE_SHLIB_HANDLE dll_handle = 0;

  ACE_TCHAR dllname[1024];
  ACE_OS::sprintf(dllname, "%s%s",ACE_DLL_PREFIX, DLL);

  ACE_TCHAR factoryname[1024];
  ACE_OS::sprintf(factoryname, "make_%s", component_name);


  dll = dllmgr->open_dll (dllname, ACE_DEFAULT_SHLIB_MODE, dll_handle );
  

  if (!dll) {
     GADGET_DEBUG1("Failed to load DLL, Possible reasons: \n");
    GADGET_DEBUG1("   * Name of DLL is wrong in XML file \n");
    GADGET_DEBUG1("   * Path of DLL is not in your DLL search path (LD_LIBRARY_PATH on Unix)\n");
    GADGET_DEBUG1("   * Path of other DLLs that this DLL depends on is not in the search path\n");
    return 0;
  }

  //Function pointer
  typedef T* (*ComponentCreator) (void);
  
  
  void *void_ptr = dll->symbol (factoryname);
  ptrdiff_t tmp = reinterpret_cast<ptrdiff_t> (void_ptr);
  ComponentCreator cc = reinterpret_cast<ComponentCreator> (tmp);

  if (cc == 0) {
    GADGET_DEBUG2("Failed to load factory (%s) from DLL (%s)\n", dllname, factoryname);
    return 0;
  }


  T* c = cc();
  
  if (!c) {
    GADGET_DEBUG1("Failed to create component using factory\n");
    return 0;
  }

  return c;
}
