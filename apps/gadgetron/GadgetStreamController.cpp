#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
#include "ace/OS_NS_stdio.h"
#include "ace/DLL.h"
#include "ace/DLL_Manager.h"
#include "ace/OS_NS_netdb.h"

#include "GadgetStreamController.h"
#include "GadgetContainerMessage.h"
#include "GadgetMessageInterface.h"
#include "GadgetronConnector.h"
#include "Gadget.h"
#include "EndGadget.h"
#include "gadgetron_config.h"

#include "gadgetron_xml.h"
#include "url_encode.h"

#include <complex>
#include <fstream>

using namespace Gadgetron;

/*
namespace Gadgetron {
//This function is needed to avoid some linking problems.
  extern "C" {
    Gadget* find_gadget_in_controller(GadgetStreamController* c, const char* g)
    {
      return c->find_gadget(g);
    }
  }
}
*/

GadgetStreamController::GadgetStreamController()
  : stream_configured_(false)
  , notifier_ (0, this, ACE_Event_Handler::WRITE_MASK)
  , writer_task_(&this->peer())
{
  gadgetron_home_ = get_gadgetron_home();
}

int GadgetStreamController::open (void)
{
	//We will set up the controllers message queue such that when a packet is enqueued write will be triggered.
	this->notifier_.reactor (this->reactor ());
	this->msg_queue ()->notification_strategy (&this->notifier_);

	ACE_TCHAR peer_name[MAXHOSTNAMELEN];
	ACE_INET_Addr peer_addr;
	if (peer().get_remote_addr (peer_addr) == 0 &&
	    peer_addr.addr_to_string (peer_name, MAXHOSTNAMELEN) == 0) {
	  GINFO("Connection from %s\n", peer_name);
	}

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
		Gadget* eg = new EndGadget();
		if (eg) {
			eg->set_controller(this);
		}
		
		ACE_NEW_RETURN(tail,
			       ACE_Module<ACE_MT_SYNCH>( ACE_TEXT("EndGadget"),
							 eg ),
			       -1);

		stream_.open(0,head,tail);
	}

	this->writer_task_.open();

	return this->reactor ()->register_handler(this,
			ACE_Event_Handler::READ_MASK);// | ACE_Event_Handler::WRITE_MASK);
}


int GadgetStreamController::handle_input (ACE_HANDLE)
{
	//Reading sequence:
	GadgetMessageIdentifier id;
	ssize_t recv_cnt = 0;
	if ((recv_cnt = peer().recv_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
	  GERROR("GadgetStreamController, unable to read message identifier\n");
	  return -1;
	}

	if (id.id == GADGET_MESSAGE_CLOSE) {
	  GDEBUG("Received close signal from client. Closing stream...\n");
	  stream_.close(1); //Shutdown gadgets and wait for them
	  GDEBUG("Stream closed\n");
	  GDEBUG("Closing writer task\n");
	  this->writer_task_.close(1);
	  GDEBUG("Writer task closed\n");
	  return 0;
	}

	GadgetMessageReader* r = readers_.find(id.id);

	if (!r) {
	  GERROR("Unrecognized Message ID received: %d\n", id.id);
	  return GADGET_FAIL;
	}

	ACE_Message_Block* mb = r->read(&peer());

	if (!mb) {
	  GERROR("GadgetMessageReader returned null pointer\n");
	  return GADGET_FAIL;
	}

	//We need to handle some special cases to make sure that we can get a stream set up.
	if (id.id == GADGET_MESSAGE_CONFIG_FILE) {
	  GadgetContainerMessage<GadgetMessageConfigurationFile>* cfgm =
	    AsContainerMessage<GadgetMessageConfigurationFile>(mb);

	  if (!cfgm) {
	    GERROR("Failed to cast message block to configuration file\n");
	    mb->release();
	    return GADGET_FAIL;
	  } else {
	    if (this->configure_from_file(std::string(cfgm->getObjectPtr()->configuration_file)) != GADGET_OK) {
	      GERROR("GadgetStream configuration failed\n");
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
	    GERROR("GadgetStream configuration failed\n");
	    mb->release();
	    return GADGET_FAIL;
	  } else {
	    mb->release();
	    return GADGET_OK;
	  }
	}

	ACE_Time_Value wait = ACE_OS::gettimeofday() + ACE_Time_Value(0,10000); //10ms from now
	if (stream_.put(mb) == -1) {
	  GERROR("Failed to put stuff on stream, too long wait, %d\n",  ACE_OS::last_error () ==  EWOULDBLOCK);
	  mb->release();
	  return GADGET_FAIL;
	}

	return GADGET_OK;
}


int GadgetStreamController::output_ready(ACE_Message_Block* mb) 
{ 
	int res = this->writer_task_.putq(mb);
	return res;
}



int GadgetStreamController::handle_close (ACE_HANDLE, ACE_Reactor_Mask mask)
{
  GDEBUG("handle_close called\n");
  
  if (mask == ACE_Event_Handler::WRITE_MASK)
    return 0;

  GINFO("Shutting down stream and closing up shop...\n");
  
  this->stream_.close();
  
  mask = ACE_Event_Handler::ALL_EVENTS_MASK |
    ACE_Event_Handler::DONT_CALL;
  
  this->reactor ()->remove_handler (this, mask);
  
  //Empty output queue in case there is something on it.
  int messages_dropped = this->msg_queue ()->flush();
  
  if (messages_dropped) {
    GDEBUG("Flushed %d messages from output queue\n", messages_dropped);
    this->reactor ()->handle_events(); //Flush any remaining events before we delete this Stream Controller
  }

  // Remove all readers and writers
  //writers_.clear();
  readers_.clear();
  
  //Clear DLL handles (to make DLLs unload if needed)
  for (unsigned int i = 0; i < dll_handles_.size(); i++) {
#if defined WIN32
    dll_handles_[i]->close(0); //On windows we will not unload the DLLs even when there are no more refs
#else 
    dll_handles_[i]->close(0); //On Unix/Mac it seems to be OK to do this
#endif
  }
  dll_handles_.clear();
  
  GINFO("Stream is closed\n");

  delete this;
  return 0;
}

Gadget* GadgetStreamController::find_gadget(std::string gadget_name)
{
  GadgetModule* gm = stream_.find(gadget_name.c_str());
  
  if (gm) {
    Gadget* g = dynamic_cast<Gadget*>(gm->writer());
		return g;
  } else {
    GDEBUG("Gadget with name %s not found! Returning null pointer\n", gadget_name.c_str());
  }

  return 0;
}

void GadgetStreamController::set_global_gadget_parameters(const std::map<std::string, std::string>& globalGadgetPara)
{
  global_gadget_parameters_ = globalGadgetPara;
}

int GadgetStreamController::configure_from_file(std::string config_xml_filename)
{
  ACE_TCHAR config_file_name[4096];
  ACE_OS::sprintf(config_file_name, "%s/%s/%s", gadgetron_home_.c_str(), GADGETRON_CONFIG_PATH, config_xml_filename.c_str());
  
  GINFO("Running configuration: %s\n", config_file_name);

  std::ifstream file (config_file_name, std::ios::in|std::ios::binary|std::ios::ate);
  if (file.is_open()) {
    size_t size = file.tellg();
    char* buffer = new char [size];
    if (!buffer) {
      GERROR("Unable to create temporary buffer for configuration file\n");
      return GADGET_FAIL;
    }
    file.seekg (0, std::ios::beg);
    file.read (buffer, size);
    file.close();
    std::string xml_file_contents(buffer,size);
    
    return configure(xml_file_contents);
    delete[] buffer;
    
  } else {
    GERROR("Unable to open configuation file: %s\n", config_file_name);
    return GADGET_FAIL;
  }
  
  return GADGET_OK;
}

int GadgetStreamController::configure(std::string config_xml_string)
{

  GadgetronXML::GadgetStreamConfiguration cfg;
  try {
    deserialize(config_xml_string.c_str(), cfg);  
  }  catch (const std::runtime_error& e) {
    GERROR("Failed to parse Gadget Stream Configuration: %s\n", e.what());
    return GADGET_FAIL;
  }

  GINFO("Found %d readers\n", cfg.reader.size());
  GINFO("Found %d writers\n", cfg.writer.size());
  GINFO("Found %d gadgets\n", cfg.gadget.size());
  
  //Configuration of readers
  for (std::vector<GadgetronXML::Reader>::iterator i = cfg.reader.begin();
       i != cfg.reader.end();
       ++i) 
    {

      long slot = 0;
      std::string dllname("");
      std::string classname("");

      slot = i->slot;
      dllname = i->dll;
      classname = i->classname;

      GINFO("--Found reader declaration\n");
      GINFO("  Reader dll: %s\n", dllname.c_str());
      GINFO("  Reader class: %s\n", classname.c_str());
      GINFO("  Reader slot: %d\n", slot);

      GadgetMessageReader* r =
	load_dll_component<GadgetMessageReader>(dllname.c_str(),
						classname.c_str());
      
      if (!r) {
	GERROR("Failed to load GadgetMessageReader from DLL\n");
	return GADGET_FAIL;
      }
      
      readers_.insert(slot, r);
      
    }	
  //Configuration of readers end


  //Configuration of writers
  for (std::vector<GadgetronXML::Writer>::iterator i = cfg.writer.begin();
       i != cfg.writer.end();
       ++i) 
    {
      long slot = 0;
      std::string dllname("");
      std::string classname("");
      
      slot = i->slot;
      dllname = i->dll;
      classname = i->classname;

      GINFO("--Found writer declaration\n");
      GINFO("  Reader dll: %s\n", dllname.c_str());
      GINFO("  Reader class: %s\n", classname.c_str());
      GINFO("  Reader slot: %d\n", slot);
      
      GadgetMessageWriter* w =
	load_dll_component<GadgetMessageWriter>(dllname.c_str(),
						classname.c_str());
      
      if (!w) {
	GERROR("Failed to load GadgetMessageWriter from DLL\n");
	return GADGET_FAIL;
      }
      
      writer_task_.register_writer(slot, w);
    }
  //Configuration of writers end

  //Let's configure the stream
  GDEBUG("Processing %d gadgets in reverse order\n",cfg.gadget.size());

  for (std::vector<GadgetronXML::Gadget>::reverse_iterator i = cfg.gadget.rbegin();
       i != cfg.gadget.rend();
       ++i) 
    {
      std::string gadgetname("");
      std::string dllname("");
      std::string classname("");

      gadgetname = i->name;
      dllname = i->dll;
      classname = i->classname;

      GINFO("--Found gadget declaration\n");
      GINFO("  Gadget Name: %s\n", gadgetname.c_str());
      GINFO("  Gadget dll: %s\n", dllname.c_str());
      GINFO("  Gadget class: %s\n", classname.c_str());

      GadgetModule* m = create_gadget_module(dllname.c_str(),
					     classname.c_str(),
					     gadgetname.c_str());
      
      if (!m) {
	GERROR("Failed to create GadgetModule from %s:%s\n",
	       classname.c_str(),
	       dllname.c_str());
	return GADGET_FAIL;
      }
      
      Gadget* g = dynamic_cast<Gadget*>(m->writer());//Get the gadget out of the module
      
      GINFO("  Gadget parameters: %d\n", i->property.size());
      for (std::vector<GadgetronXML::GadgetronParameter>::iterator p = i->property.begin();
	   p != i->property.end();
	   ++p)
	{
	  std::string pname(p->name);
	  std::string pval(p->value);
	  GINFO("Setting parameter %s = %s\n", pname.c_str(),pval.c_str());
	  g->set_parameter(pname.c_str(),pval.c_str(),false);
	}
      
        // set the global gadget parameters for every gadget
      std::map<std::string, std::string>::const_iterator iter;
      for ( iter=global_gadget_parameters_.begin(); iter!=global_gadget_parameters_.end(); iter++ )
        {
	  std::string key = iter->first;
	  std::string value = iter->second;
	  g->set_parameter(key.c_str(), value.c_str(), false);
        }

      if (stream_.push(m) < 0) {
	GERROR("Failed to push Gadget %s onto stream\n", gadgetname.c_str());
	delete m;
	return GADGET_FAIL;
      }
      
    }

  GINFO("Gadget Stream configured\n");
  stream_configured_ = true;

  return GADGET_OK;
}

GadgetModule * GadgetStreamController::create_gadget_module(const char* DLL, 
		const char* gadget,
		const char* gadget_module_name)
{

	Gadget* g = load_dll_component<Gadget>(DLL,gadget);

	if (!g) {
	  GERROR("Failed to load gadget using factory\n");
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
#if defined(WIN32) && defined(_DEBUG)
	ACE_OS::sprintf(dllname, "%s%sd",ACE_DLL_PREFIX, DLL);
#else
	ACE_OS::sprintf(dllname, "%s%s",ACE_DLL_PREFIX, DLL);
#endif

	ACE_TCHAR factoryname[1024];
	ACE_OS::sprintf(factoryname, "make_%s", component_name);

	dll = dllmgr->open_dll (dllname, ACE_DEFAULT_SHLIB_MODE, dll_handle );

	if (!dll) {
	  GERROR("Failed to load DLL, Possible reasons: \n");
	  GERROR("   * Name of DLL is wrong in XML file \n");
	  GERROR("   * Path of DLL is not in your DLL search path (LD_LIBRARY_PATH on Unix)\n");
	  GERROR("   * Path of other DLLs that this DLL depends on is not in the search path\n");
	  return 0;
	} else {
	  dll_handles_.push_back(dll);
	}

	//Function pointer
	typedef T* (*ComponentCreator) (void);

	void *void_ptr = dll->symbol (factoryname);
	ptrdiff_t tmp = reinterpret_cast<ptrdiff_t> (void_ptr);
	ComponentCreator cc = reinterpret_cast<ComponentCreator> (tmp);

	if (cc == 0) {
	  GERROR("Failed to load factory (%s) from DLL (%s)\n", dllname, factoryname);
	  return 0;
	}

	T* c = cc();

	if (!c) {
	  GERROR("Failed to create component using factory\n");
	  return 0;
	}

	return c;
}
