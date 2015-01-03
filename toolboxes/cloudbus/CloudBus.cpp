#include "CloudBus.h"
#include "log.h"

namespace Gadgetron
{
  CloudBus* CloudBus::instance_ = 0;
  const char* CloudBus::mcast_inet_addr_ = GADGETRON_DEFAULT_MULTICAST_ADDR;
  int CloudBus::mcast_port_ = GADGETRON_DEFAULT_MULTICAST_PORT;
  bool CloudBus::query_mode_ = false; //Listen only is disabled default
  int CloudBus::gadgetron_port_ = 9002; //Default port

  CloudBusTask::CloudBusTask(int port, const char* addr)
    : inherited()
    , mcast_addr_(port, addr)
    , mcast_dgram_(ACE_SOCK_Dgram_Mcast::OPT_BINDADDR_NO)
  {
  }

  CloudBusTask::CloudBusTask()
    : inherited()
    , mcast_addr_(GADGETRON_DEFAULT_MULTICAST_PORT, GADGETRON_DEFAULT_MULTICAST_ADDR)
    , mcast_dgram_(ACE_SOCK_Dgram_Mcast::OPT_BINDADDR_NO)
  {
  }
    
  int CloudBusTask::open(void*)
  {
    return this->activate( THR_NEW_LWP | THR_JOINABLE,1); //single thread
  }

  CloudBusReceiverTask::CloudBusReceiverTask(int port, const char* addr)
    : CloudBusTask(port, addr)
  {
    
  }

  int CloudBusReceiverTask::open(void*)
  {
    if (mcast_dgram_.join(mcast_addr_) == -1) {
      GDEBUG_STREAM("Error doing dgram join" << std::endl);
      return -1;
    }
    return CloudBusTask::open();      
  }

  int CloudBusReceiverTask::close(u_long flags)
  {
    mcast_dgram_.leave(mcast_addr_);
    return CloudBusTask::close(flags);
  }

  int CloudBusReceiverTask::svc(void)
  {
    char buffer[GADGETRON_NODE_INFO_MESSAGE_LENGTH]; //Size of message
    GadgetronNodeInfo info;
    ACE_INET_Addr peer_address;
    while (mcast_dgram_.recv(buffer, GADGETRON_NODE_INFO_MESSAGE_LENGTH, peer_address) != -1)
      {
	info.uuid = boost::uuids::to_string(*((boost::uuids::uuid*)buffer));
	info.address = std::string(peer_address.get_host_addr());
	memcpy(&info.port              , buffer + 16,                    sizeof(uint32_t));
	memcpy(&info.compute_capability, buffer + 16 + sizeof(uint32_t), sizeof(uint32_t));
	CloudBus::instance()->update_node(info.uuid.c_str(), info);
      }

    return 0;
  }

  CloudBusSenderTask::CloudBusSenderTask(int port, const char* addr)
    : CloudBusTask(port, addr)
  {
    
  }

  int CloudBusSenderTask::open(void*)
  {
    if (mcast_dgram_.open(mcast_addr_) == -1) {
      GDEBUG_STREAM("Error doing dgram open" << std::endl);
      return -1;
    }
    return CloudBusTask::open();      
  }

  int CloudBusSenderTask::svc(void)
  {
    char buffer[GADGETRON_NODE_INFO_MESSAGE_LENGTH]; //Size of message
    if (CloudBus::instance()->uuid_.size() != 16) {
      GDEBUG_STREAM("Severe problem, UUID is != 16" << std::endl);
      GDEBUG_STREAM("uuid: " << CloudBus::instance()->uuid_ << "(" << CloudBus::instance()->uuid_.size() << ")" << std::endl);
    }
    
    memcpy(buffer                        ,  CloudBus::instance()->uuid_.begin(), 16);
    memcpy(buffer + 16,                    &CloudBus::instance()->node_info_.port, sizeof(uint32_t));
    memcpy(buffer + 16 + sizeof(uint32_t), &CloudBus::instance()->node_info_.compute_capability, sizeof(uint32_t));
    
    while (true) {
      if (!CloudBus::instance()->query_mode_) {
	if (mcast_dgram_.send(buffer, GADGETRON_NODE_INFO_MESSAGE_LENGTH) == -1) {
	  GDEBUG_STREAM("Failed to send dgram data" << std::endl);
	}
      }
      CloudBus::instance()->remove_stale_nodes();
      ACE_OS::sleep(5);//Sleep for 5 seconds
    }

    return 0;
  }

  CloudBus* CloudBus::instance()
  {
    if (!instance_)
      {
	instance_ = new CloudBus(mcast_port_, mcast_inet_addr_);
	instance_->receiver_.open();
	instance_->sender_.open();
      }
    return instance_;
  }

  
  void CloudBus::set_mcast_address(const char* addr)
  {
    mcast_inet_addr_ = addr;
  }

  void CloudBus::set_mcast_port(int port)
  {
    mcast_port_ = port;
  }

  void CloudBus::set_query_only(bool m)
  {
    query_mode_ = m;
  }
 
  void CloudBus::set_gadgetron_port(uint32_t port)
  {
    gadgetron_port_ = port;
  }

  void CloudBus::wait()
  {
    sender_.wait();
    receiver_.wait();
    receiver_.close();
  }

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

  CloudBus::CloudBus(int port, const char* addr)
    : receiver_(port, addr)
    , sender_(port, addr)
    , mtx_("CLOUDBUSMTX")
    , uuid_(boost::uuids::random_generator()())
  {
    node_info_.port = gadgetron_port_;
    set_compute_capability(1);
    node_info_.uuid = boost::uuids::to_string(uuid_);
    ACE_SOCK_Acceptor listener (ACE_Addr::sap_any);
    ACE_INET_Addr local_addr;
    listener.get_local_addr (local_addr);
    node_info_.address = std::string(local_addr.get_host_name());
  }

  void CloudBus::update_node(const char* a, GadgetronNodeInfo& info)
  {
    mtx_.acquire();
    std::string key(a);
    map_type_::iterator it = nodes_.find(key);
    if (it == nodes_.end()) {
      if (info.uuid != node_info_.uuid) { //Reject stuff coming from myself
	GDEBUG_STREAM("---->>>> New Cloud Node <<<<< ----- " << info.uuid << " (" << info.address << ":" << info.port << ", " << info.compute_capability << ")" << std::endl);
      } 
    } 

    if (info.uuid != node_info_.uuid) {
      nodes_[key] = std::pair<GadgetronNodeInfo,time_t>(info,time(NULL));
    }
    mtx_.release();
  }

  void CloudBus::remove_stale_nodes()
  {
    mtx_.acquire();
    map_type_ new_nodes_;
    time_t now = time(NULL);
    for (map_type_::iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
      if (fabs(difftime(it->second.second,now)) > 30) {
        GadgetronNodeInfo n = it->second.first;
        GDEBUG_STREAM("---->>>> DELETING STALE CLOUD NODE <<<<< ----- " << n.uuid << " (" << n.address << ":" << n.port  << ", " << n.compute_capability << ")" << std::endl);
      }
      else
      {
        new_nodes_[it->first] = it->second;
      }
    }

    nodes_.clear();
    nodes_ = new_nodes_;

    mtx_.release();
  }
}
