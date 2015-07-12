#include "GadgetronConnector.h"

#include <ace/SOCK_Connector.h>
#include "log.h"

using namespace Gadgetron;

GadgetronConnector::GadgetronConnector()
  : writer_task_(&this->peer())
{
}

GadgetronConnector::~GadgetronConnector() {
  readers_.clear();
}

int GadgetronConnector::open(std::string hostname, std::string port)
{
  hostname_= hostname;
  port_ = port;

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

  this->writer_task_.open();

  return this->activate( THR_NEW_LWP | THR_JOINABLE, 1);
}

int GadgetronConnector::handle_input(ACE_HANDLE fd)
{
  return 0;
}

int GadgetronConnector::handle_close(ACE_HANDLE handle, ACE_Reactor_Mask mask)
{

  GDEBUG(" >>>>>>>>>>>>> HANDLE_CLOSE <<<<<<<<<<<<<<<<\n");
  if (mask == ACE_Event_Handler::WRITE_MASK)
    return 0;

  peer().close();

  return 0;
}

int GadgetronConnector::svc(void)
{
  while (true) {
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
  }
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

  if (this->peer().send_n(config_xml.c_str(), ini.script_length) != ini.script_length) {
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

  if (this->peer().send_n(xml_string.c_str(), conf.script_length) != conf.script_length) {
    GERROR("Unable to send parameter xml\n");
    return -1;
  }
  return 0;
}

