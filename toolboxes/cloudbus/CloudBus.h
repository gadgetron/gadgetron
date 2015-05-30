#ifndef GADGETRON_CLOUDBUS_H
#define GADGETRON_CLOUDBUS_H

#include "cloudbus_export.h"
#include <ace/Task.h>
#include <ace/INET_Addr.h>
#include <ace/OS_NS_unistd.h>
#include <ace/SOCK_Connector.h>
#include <ace/SOCK_Stream.h>
#include <ace/SOCK_Acceptor.h>
#include <ace/Svc_Handler.h>

#include "log.h"

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <iostream>
#include <map>
#include <utility>
#include <time.h>
#include <vector>

#define GADGETRON_DEFAULT_RELAY_ADDR "localhost"
#define GADGETRON_DEFAULT_RELAY_PORT 8002
#define GADGETRON_NODE_INFO_MESSAGE_LENGTH 16+sizeof(uint32_t)*2 //16 bytes for uuid + 2 ints

#define MAXHOSTNAMELENGTH 1024

namespace Gadgetron
{

  struct GadgetronNodeInfo
  {
    std::string uuid;
    std::string address;
    uint32_t port;
    uint32_t compute_capability;
  };
  

  class EXPORTCLOUDBUS CloudBus : public ACE_Svc_Handler<ACE_SOCK_STREAM, ACE_MT_SYNCH>
  {
  public:
    typedef ACE_Task<ACE_MT_SYNCH> inherited;    
    virtual int open(void* = 0)
    {
      if (!this->reactor()) {
	this->reactor(ACE_Reactor::instance());
      }
      return this->activate( THR_NEW_LWP | THR_JOINABLE,1); //single thread
    }

    virtual int handle_close (ACE_HANDLE handle, ACE_Reactor_Mask close_mask)
    {
      GDEBUG("Cloud bus connection closed\n");
      this->peer().close();
      connected_ = false;
      return 0;
    }

    virtual int svc(void)
    {
      
      while (true) {
	if (connected_) {
	  //
	} else {
	  std::string connect_addr(relay_inet_addr_);
	  if (connect_addr == "localhost") {
	    connect_addr = node_info_.address;
	  }
	  ACE_INET_Addr server(relay_port_,connect_addr.c_str());
	  ACE_SOCK_Connector connector;

	  if (connector.connect(this->peer(),server) == 0) {
	    ACE_TCHAR peer_name[MAXHOSTNAMELENGTH];
	    ACE_INET_Addr peer_addr;
	    if ((this->peer().get_remote_addr (peer_addr) == 0) && 
		(peer_addr.addr_to_string (peer_name, MAXHOSTNAMELENGTH) == 0)) {

	      GDEBUG("Connected to %s\n", peer_name);
	      if (this->reactor ()->register_handler(this, ACE_Event_Handler::READ_MASK) != 0) {
		GERROR("Failed to register read handler\n");
		return -1;
	      }
	      connected_ = true;
	    } 
	  }
	}
	//Sleep for 5 seconds
	ACE_Time_Value tv (5);
	ACE_OS::sleep (tv);	  	
      }
      return 0;
    }

    static CloudBus* instance();
    static void set_relay_address(const char* addr);
    static void set_relay_port(int port);
    static void set_query_only(bool m = true);
    static void set_gadgetron_port(uint32_t port);

    void set_compute_capability(uint32_t c)
    {
      node_info_.compute_capability = c;
    }


    //void get_node_info(std::vector<GadgetronNodeInfo>& nodes);
    //size_t get_number_of_nodes();

  protected:
    ///Protected constructor. 
    CloudBus(int port, const char* addr);
    CloudBus(); 

    static CloudBus* instance_;
    static const char* relay_inet_addr_;
    static int relay_port_;
    static bool query_mode_; //Listen only
    static int gadgetron_port_;

    GadgetronNodeInfo node_info_;
    
    ACE_Thread_Mutex mtx_;

    boost::uuids::uuid uuid_;
    bool connected_;
    ACE_SOCK_Stream socket_;
  };


}

#endif
