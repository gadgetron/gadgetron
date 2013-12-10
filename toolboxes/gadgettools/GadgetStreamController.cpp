#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
#include "ace/OS_NS_stdio.h"
#include "ace/DLL.h"
#include "ace/DLL_Manager.h"
#include "ace/OS_NS_netdb.h"

#include "GadgetStreamController.h"
#include "GadgetContainerMessage.h"
#include "Gadget.h"
#include "EndGadget.h"

#include "gadgetron.hxx" //Auto generated class representation of gadgetron XML configuration
#include "url_encode.h"

#include <complex>
#include <fstream>

using namespace Gadgetron;
int GadgetStreamController::open (void)
{
	//We will set up the controllers message queue such that when a packet is enqueued write will be triggered.
	this->notifier_.reactor (this->reactor ());
	this->msg_queue ()->notification_strategy (&this->notifier_);
    this->msg_queue()->high_water_mark((size_t)(48.0*1024*1024*1024));

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
		ACE_DEBUG ((LM_DEBUG,
				ACE_TEXT ("(%P|%t) GadgetStreamController, unable to read message identifier\n")));
		return -1;
	}

	if (id.id == GADGET_MESSAGE_CLOSE) {
		GADGET_DEBUG1("Received close signal from client. Closing stream...\n");
		stream_.close(1); //Shutdown gadgets and wait for them
		GADGET_DEBUG1("Stream closed\n");
		GADGET_DEBUG1("Closing writer task\n");
		this->writer_task_.close(1);
		GADGET_DEBUG1("Writer task closed\n");
		return 0;
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


int GadgetStreamController::output_ready(ACE_Message_Block* mb) 
{ 
	int res = this->writer_task_.putq(mb);
	return res;
}



int GadgetStreamController::handle_close (ACE_HANDLE, ACE_Reactor_Mask mask)
{
	GADGET_DEBUG1("handle_close called\n");

	if (mask == ACE_Event_Handler::WRITE_MASK)
		return 0;

	GADGET_DEBUG1("Shutting down stream and closing up shop...\n");

	this->stream_.close();

	mask = ACE_Event_Handler::ALL_EVENTS_MASK |
			ACE_Event_Handler::DONT_CALL;

	this->reactor ()->remove_handler (this, mask);

	//Empty output queue in case there is something on it.
	int messages_dropped = this->msg_queue ()->flush();

	if (messages_dropped) {
		GADGET_DEBUG2("Flushed %d messages from output queue\n", messages_dropped);
		this->reactor ()->handle_events(); //Flush any remaining events before we delete this Stream Controller
	}

	// Remove all readers and writers
	//writers_.clear();
	readers_.clear();

	//Clear DLL handles (to make DLLs unload if needed)
	for (size_t i = 0; i < dll_handles_.size(); i++) {
#if defined WIN32
		dll_handles_[i]->close(0); //On windows we will not unload the DLLs even when there are no more refs
#else 
		dll_handles_[i]->close(0); //On Unix/Mac it seems to be OK to do this
#endif
	}
	dll_handles_.clear();

	GADGET_DEBUG1("Stream is closed\n");

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
		GADGET_DEBUG2("Gadget with name %s not found! Returning null pointer\n", gadget_name.c_str());
	}

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

		return configure(xml_file_contents);
		delete[] buffer;

	} else {
		GADGET_DEBUG2("Unable to open configuation file: %s\n", config_file_name);
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

int GadgetStreamController::configure(std::string config_xml_string)
{

	char * gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");
	ACE_TCHAR schema_file_name[4096];
	ACE_OS::sprintf(schema_file_name, "%s/schema/gadgetron.xsd", gadgetron_home);

	std::string tmp(schema_file_name);
	tmp = url_encode(tmp);
	ACE_OS_String::strncpy(schema_file_name,tmp.c_str(), 4096);


	xml_schema::properties props;
	props.schema_location (
	  "http://gadgetron.sf.net/gadgetron",
	  std::string (schema_file_name));

	std::istringstream str_stream(config_xml_string, std::stringstream::in);
	std::auto_ptr<gadgetron::gadgetronStreamConfiguration> cfg;

	try {
		cfg = std::auto_ptr<gadgetron::gadgetronStreamConfiguration>(gadgetron::gadgetronStreamConfiguration_(str_stream,0,props));
		//cfg = std::auto_ptr<gadgetron::gadgetronStreamConfiguration>(gadgetron::gadgetronStreamConfiguration_(std::string(config_file_name)));
	}  catch (const xml_schema::exception& e) {
		GADGET_DEBUG2("Failed to parse Gadget Stream Configuration: %s\n", e.what());
		return GADGET_FAIL;
	}

	GADGET_DEBUG2("Found %d readers\n", cfg->reader().size());
	GADGET_DEBUG2("Found %d writers\n", cfg->writer().size());
	GADGET_DEBUG2("Found %d gadgets\n", cfg->gadget().size());

	for (gadgetron::gadgetronStreamConfiguration::reader_sequence::iterator i (cfg->reader().begin ()); i != cfg->reader().end(); ++i) {
		long slot = 0;
		std::string dllname("");
		std::string classname("");

		slot = i->slot();
		dllname = i->dll();
		classname = i->classname();

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

		readers_.insert((unsigned short)slot, r);

	}
	//Configuration of readers end


	//Configuration of writers
	for (gadgetron::gadgetronStreamConfiguration::writer_sequence::iterator i (cfg->writer().begin ()); i != cfg->writer().end(); ++i) {
		long slot = 0;
		std::string dllname("");
		std::string classname("");

		slot = i->slot();
		dllname = i->dll();
		classname = i->classname();

		GADGET_DEBUG1("--Found writer declaration\n");
		GADGET_DEBUG2("  Reader dll: %s\n", dllname.c_str());
		GADGET_DEBUG2("  Reader class: %s\n", classname.c_str());
		GADGET_DEBUG2("  Reader slot: %d\n", slot);

		GadgetMessageWriter* w =
				load_dll_component<GadgetMessageWriter>(dllname.c_str(),
						classname.c_str());

		if (!w) {
			GADGET_DEBUG1("Failed to load GadgetMessageWriter from DLL\n");
			return GADGET_FAIL;
		}

		writer_task_.register_writer(slot, w);
	}
	//Configuration of writers end

	//Let's configure the stream
	GADGET_DEBUG2("Processing %d gadgets in reverse order\n",cfg->gadget().size());
	for (gadgetron::gadgetronStreamConfiguration::gadget_sequence::reverse_iterator i (cfg->gadget().rbegin ()); i != cfg->gadget().rend(); ++i) {
		std::string gadgetname("");
		std::string dllname("");
		std::string classname("");

		gadgetname = i->name();
		dllname = i->dll();
		classname = i->classname();

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

		Gadget* g = dynamic_cast<Gadget*>(m->writer());//Get the gadget out of the module

		GADGET_DEBUG2("  Gadget parameters: %d\n", i->property().size());
		for (gadgetron::gadget::property_sequence::iterator p (i->property().begin()); p != i->property().end(); ++p) {
			std::string pname(p->name());
			std::string pval(p->value());
			GADGET_DEBUG2("Setting parameter %s = %s\n", pname.c_str(),pval.c_str());
			g->set_parameter(pname.c_str(),pval.c_str(),false);
		}

		if (stream_.push(m) < 0) {
			GADGET_DEBUG2("Failed to push Gadget %s onto stream\n", gadgetname.c_str());
			delete m;
			return GADGET_FAIL;
		}

	}

	GADGET_DEBUG1("Gadget Stream configured\n");
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
#if defined(WIN32) && defined(_DEBUG)
	ACE_OS::sprintf(dllname, "%s%sd",ACE_DLL_PREFIX, DLL);
#else
	ACE_OS::sprintf(dllname, "%s%s",ACE_DLL_PREFIX, DLL);
#endif

	ACE_TCHAR factoryname[1024];
	ACE_OS::sprintf(factoryname, "make_%s", component_name);

	dll = dllmgr->open_dll (dllname, ACE_DEFAULT_SHLIB_MODE, dll_handle );

	if (!dll) {
		GADGET_DEBUG1("Failed to load DLL, Possible reasons: \n");
		GADGET_DEBUG1("   * Name of DLL is wrong in XML file \n");
		GADGET_DEBUG1("   * Path of DLL is not in your DLL search path (LD_LIBRARY_PATH on Unix)\n");
		GADGET_DEBUG1("   * Path of other DLLs that this DLL depends on is not in the search path\n");
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
