#include "CloudBus.h"
#include "log.h"

namespace Gadgetron
{
  CloudBus* CloudBus::instance_ = 0;
  const char* CloudBus::relay_inet_addr_ = GADGETRON_DEFAULT_RELAY_ADDR;
  int CloudBus::relay_port_ = GADGETRON_DEFAULT_RELAY_PORT;
  bool CloudBus::query_mode_ = false; //Listen only is disabled default
  int CloudBus::gadgetron_port_ = 9002; //Default port

  CloudBus* CloudBus::instance()
  {
    if (!instance_)
      {
	instance_ = new CloudBus(relay_port_, relay_inet_addr_);
	instance_->open();
      }
    return instance_;
  }

  
  void CloudBus::set_relay_address(const char* addr)
  {
    relay_inet_addr_ = addr;
  }

  void CloudBus::set_relay_port(int port)
  {
    relay_port_ = port;
  }

  void CloudBus::set_query_only(bool m)
  {
    query_mode_ = m;
  }
 
  void CloudBus::set_gadgetron_port(uint32_t port)
  {
    gadgetron_port_ = port;
  }

  /*
  void CloudBus::get_node_info(std::vector<GadgetronNodeInfo>& nodes)
  {
    mtx_.acquire();
    nodes.clear();
    for (map_type_::iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
      GadgetronNodeInfo n = it->second.first;
      nodes.push_back(n);
    }
    mtx_.release();
  }
  
  size_t CloudBus::get_number_of_nodes()
  {
    size_t n = 0;
    mtx_.acquire();
    n = nodes_.size();
    mtx_.release();
    return n;
  }
  */

  CloudBus::CloudBus(int port, const char* addr)
    : mtx_("CLOUDBUSMTX")
    , uuid_(boost::uuids::random_generator()())
    , connected_(false)
  {
    node_info_.port = gadgetron_port_;
    set_compute_capability(1);
    node_info_.uuid = boost::uuids::to_string(uuid_);
    ACE_SOCK_Acceptor listener (ACE_Addr::sap_any);
    ACE_INET_Addr local_addr;
    listener.get_local_addr (local_addr);
    node_info_.address = std::string(local_addr.get_host_name());
  }

}
