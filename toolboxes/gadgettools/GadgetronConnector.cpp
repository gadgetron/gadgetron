#include "GadgetronConnector.h"

#include <ace/SOCK_Connector.h>
#include "log.h"

using namespace Gadgetron;

GadgetronConnector::GadgetronConnector()
    //: notifier_ (0, this, ACE_Event_Handler::WRITE_MASK)
    : writer_task_(&this->peer())
{
}

GadgetronConnector::~GadgetronConnector() {
    readers_.clear();
    //writers_.clear();
}

int GadgetronConnector::openImpl(std::string hostname, std::string port)
{
    hostname_= hostname;
    port_ = port;

    //We will add a notification strategy to the message queue to make sure than handle_output gets triggered when packages are on the queue
    //this->notifier_.reactor (this->reactor ());
    //this->msg_queue ()->notification_strategy (&this->notifier_);

    ACE_INET_Addr server(port_.c_str(),hostname_.c_str());
    ACE_SOCK_Connector connector;

    if (connector.connect(this->peer(),server) == -1) {
      GERROR("Failed to connect");
      return -1;
    }

    ACE_TCHAR peer_name[MAXHOSTNAMELENGTH];
    ACE_INET_Addr peer_addr;
    if (peer().get_remote_addr (peer_addr) == 0 && peer_addr.addr_to_string (peer_name, MAXHOSTNAMELENGTH) == 0) {
      GDEBUG("Connection from %s\n", peer_name);
    }

    return 0;
}

int GadgetronConnector::open(std::string hostname, std::string port)
{
    //Make sure we have a reactor, otherwise assign one from the singleton instance
    if (!this->reactor()) {
      GDEBUG("Setting reactor");
      this->reactor(ACE_Reactor::instance());
    }

    this->openImpl(hostname, port);

    this->writer_task_.open();

    if (this->reactor ()->register_handler(this, ACE_Event_Handler::READ_MASK) != 0) {
      GERROR("Failed to register read handler\n");
      return -2;
    }

    return this->activate( THR_NEW_LWP | THR_JOINABLE, 1); //Run single threaded. TODO: Add multithreaded support
}

int GadgetronConnector::handle_input(ACE_HANDLE /*fd*/)
{
    ssize_t recv_count = 0;
    GadgetMessageIdentifier mid;

    if ((recv_count = peer().recv_n(&mid, sizeof(GadgetMessageIdentifier))) <= 0) {
      GWARN("GadgetronConnector, failed to read message identifier\n");
      return -1;
    }

    //Is this a shutdown message?
    if (mid.id == GADGET_MESSAGE_CLOSE) {
      GDEBUG("GadgetronConnector, Close Message received\n");
      return close();
    }

    GadgetMessageReader* r = readers_.find(mid.id);
    if (r == 0) {
      GERROR("GadgetronConnector, Unknown message id %d received\n", mid.id);
      return -1;
    }

    ACE_Message_Block* mb = r->read(&peer());

    if (!mb) {
      GERROR("GadgetronConnector, Failed to read message\n");
      return -1;
    }  else {
      if (process(mid.id, mb) < 0) {
	GERROR("GadgetronConnector, Failed to process message\n");
	return -1;
      }
    }

    return 0;
}

int GadgetronConnector::handle_close(ACE_HANDLE /*handle*/, ACE_Reactor_Mask /*close_mask*/)
{
  GDEBUG("Handling close...\n");
  this->reactor()->end_reactor_event_loop();
  return 0;//this->wait();
}

int GadgetronConnector::svc(void)
{
    //ACE_thread_t old_owner;

    //Take ownership of Reactor
    this->reactor()->owner(ACE_Thread::self ());//, &old_owner);

    this->reactor()->reset_event_loop();

    ACE_Time_Value initialDelay (3);
    ACE_Time_Value interval (0,100);

    //Handle the events
    this->reactor()->run_reactor_event_loop();

    //this->reactor()->owner(&old_owner);
    
    GDEBUG("GadgetronConnector svc done...\n");

    return 0;
}

int GadgetronConnector::register_reader(size_t slot, GadgetMessageReader *reader)
{
    return readers_.insert( (unsigned short)slot,reader);
}

int GadgetronConnector::send_gadgetron_configuration_file(std::string config_xml_name)
{
    GadgetMessageIdentifier id;
    id.id = GADGET_MESSAGE_CONFIG_FILE;

    GadgetMessageConfigurationFile ini;
    ACE_OS_String::strncpy(ini.configuration_file, config_xml_name.c_str(),1024);


    if (this->peer().send_n(&id, sizeof(GadgetMessageIdentifier)) != sizeof(GadgetMessageIdentifier)) {
      GERROR("Unable to send GadgetMessageIdentifier\n");
      return -1;
    }

    if (this->peer().send_n(&ini, sizeof(GadgetMessageConfigurationFile)) != sizeof(GadgetMessageConfigurationFile)) {
      GERROR("Unable to send GadgetMessageConfigurationFile\n");
      return -1;
    }

    return 0;
}

int GadgetronConnector::send_gadgetron_configuration_script(std::string config_xml)
{
    GadgetMessageIdentifier id;
    id.id = GADGET_MESSAGE_CONFIG_SCRIPT;

    GadgetMessageScript ini;
    ini.script_length = (ACE_UINT32)config_xml.size()+1;

    if (this->peer().send_n(&id, sizeof(GadgetMessageIdentifier)) != sizeof(GadgetMessageIdentifier)) {
      GERROR("Unable to send GadgetMessageIdentifier\n");
      return -1;
    }

    if (this->peer().send_n(&ini, sizeof(GadgetMessageScript)) != sizeof(GadgetMessageScript)) {
      GERROR("Unable to send GadgetMessageScript\n");
      return -1;
    }

    if (this->peer().send_n(config_xml.c_str(), ini.script_length) != static_cast<ssize_t>(ini.script_length)) {
      GERROR("Unable to send parameter xml\n");
      return -1;
    }

    return 0;
}

int GadgetronConnector::send_gadgetron_parameters(std::string xml_string)
{
    GadgetMessageIdentifier id;
    id.id = GADGET_MESSAGE_PARAMETER_SCRIPT;

    GadgetMessageScript conf;
    conf.script_length = (ACE_UINT32)xml_string.size()+1;
    if (this->peer().send_n(&id, sizeof(GadgetMessageIdentifier)) != sizeof(GadgetMessageIdentifier)) {
      GERROR("Unable to send GadgetMessageIdentifier\n");
      return -1;
    }

    if (this->peer().send_n(&conf, sizeof(GadgetMessageScript)) != sizeof(GadgetMessageScript)) {
      GERROR("Unable to send GadgetMessageScript\n");
      return -1;
    }

    if (this->peer().send_n(xml_string.c_str(), conf.script_length) != static_cast<ssize_t>(conf.script_length)) {
      GERROR("Unable to send parameter xml\n");
      return -1;
    }

    return 0;
}

