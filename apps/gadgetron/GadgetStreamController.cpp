#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
#include "ace/OS_NS_stdio.h"
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
#include "CloudBus.h"

#include <complex>
#include <fstream>
#include <boost/filesystem.hpp>

#ifdef USE_GTBABYLON
    #include <GTBabylon.h>
#endif

using namespace Gadgetron;

GadgetStreamController::GadgetStreamController()
  : GadgetStreamInterface()
  , notifier_ (0, this, ACE_Event_Handler::WRITE_MASK)
  , writer_task_(&this->peer())
{
  CloudBus::instance()->report_recon_start();    
}

GadgetStreamController::~GadgetStreamController()
{ 
  CloudBus::instance()->report_recon_end();
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


  this->writer_task_.open();

  return this->activate( THR_NEW_LWP | THR_DETACHED, 1);

}

int GadgetStreamController::svc(void)
{
  while (true) {
    GadgetMessageIdentifier id;
    ssize_t recv_cnt = 0;
    if ((recv_cnt = peer().recv_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
      GERROR("GadgetStreamController, unable to read message identifier\n");
      return -1;
    }

    if (id.id == GADGET_MESSAGE_CLOSE) {
      stream_.close(1); //Shutdown gadgets and wait for them
      GDEBUG("Stream closed\n");
      GDEBUG("Closing writer task\n");
      this->writer_task_.close(1);
      GDEBUG("Writer task closed\n");
      continue;
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
	  continue;
	}
      }
    } else if (id.id == GADGET_MESSAGE_CONFIG_SCRIPT) {
      std::string xml_config(mb->rd_ptr(), mb->length());
      //std::stringstream stream(xml_config, std::ios::in);
      if (this->configure(xml_config) != GADGET_OK) {
	GERROR("GadgetStream configuration failed\n");
	mb->release();
	return GADGET_FAIL;
      } else {
	mb->release();
	continue;
      }
    }

    ACE_Time_Value wait = ACE_OS::gettimeofday() + ACE_Time_Value(0,10000); //10ms from now
    if (stream_.put(mb) == -1) {
      GERROR("Failed to put stuff on stream, too long wait, %d\n",  ACE_OS::last_error () ==  EWOULDBLOCK);
      mb->release();
      return GADGET_FAIL;
    }
  }
  return GADGET_OK;
}

int GadgetStreamController::handle_input (ACE_HANDLE)
{
  return 0;
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

namespace Babylon
{
#ifdef USE_GTBABYLON
    std::string verify_signature(const std::string& config_string)
    {
        //GDEBUG_STREAM(config_string);
        return GTBabylon::decode_message(config_string);
    }
#else
    std::string verify_signature(const std::string& config_string)
    {
        //GDEBUG_STREAM(config_string);
        // Config file signature verification is not enabled for this build - just pass on the config, it's fine.
        return config_string;
    }
#endif
}

int GadgetStreamController::configure_from_file(std::string filename)
{
  boost::filesystem::path full_path = gadgetron_home_ / GADGETRON_CONFIG_PATH / filename;

  GINFO("Running configuration: %s\n", full_path.c_str());

  std::ifstream config_file_stream (full_path.c_str(), std::ios::in);

  if (!config_file_stream.is_open()) {
    GERROR("Unable to open configuation file: %s\n", full_path.c_str());
    return GADGET_FAIL;
  }

  auto config_string = std::string(std::istreambuf_iterator<char>(config_file_stream), {});
  std::string config_string_used = Babylon::verify_signature(config_string);
  return configure(config_string_used);
}

int GadgetStreamController::configure(const std::string& config_string)
{
  std::string config_string_used = config_string;

  if(config_string.find("GTBABYLON")!=std::string::npos)
  {
      config_string_used = Babylon::verify_signature(config_string);
  }

  GadgetronXML::GadgetStreamConfiguration cfg;
  try {
    deserialize(config_string_used, cfg);
  }  catch (const std::runtime_error& e) {
    GERROR("Failed to parse Gadget Stream Configuration: %s\n", e.what());
    return GADGET_FAIL;
  }

  stream_configuration_ = cfg;

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
      GINFO("  Writer dll: %s\n", dllname.c_str());
      GINFO("  Writer class: %s\n", classname.c_str());
      GINFO("  Writer slot: %d\n", slot);
      
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


