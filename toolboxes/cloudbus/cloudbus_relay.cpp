#include "gadgetron_rest.h"

#include "ace/Reactor.h"
#include "ace/SOCK_Acceptor.h"
#include "ace/Svc_Handler.h"
#include "ace/SOCK_Stream.h"
#include "ace/Reactor_Notification_Strategy.h"
#include "ace/Stream.h"

#include <map>
#include <chrono>

#include "log.h"
#include "CloudBus.h"

namespace Gadgetron{

  //Forward declaration
  class CloudBusNodeController;

  class CloudBusRelayAcceptor : public ACE_Event_Handler
  {
  public:
    CloudBusRelayAcceptor()
      : mtx_("CLOUDBUSRELAYMTX")
    {

    }

    virtual ~CloudBusRelayAcceptor ()
    {
      this->handle_close (ACE_INVALID_HANDLE, 0);
    }

    int open (const ACE_INET_Addr &listen_addr)
    {
      if (this->acceptor_.open (listen_addr, 1) == -1) {
	GERROR("error opening acceptor\n");
	return -1;
      }
      return this->reactor ()->register_handler(this, ACE_Event_Handler::ACCEPT_MASK);
    }

    virtual ACE_HANDLE get_handle (void) const
    {
      return this->acceptor_.get_handle ();
    }

    virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);

    virtual int handle_close (ACE_HANDLE handle,
			      ACE_Reactor_Mask close_mask)
    {
      GDEBUG("Close CloudBus Relay Acceptor\n");

      if (this->acceptor_.get_handle () != ACE_INVALID_HANDLE) {
	ACE_Reactor_Mask m =
	  ACE_Event_Handler::ACCEPT_MASK | ACE_Event_Handler::DONT_CALL;
	this->reactor ()->remove_handler (this, m);
	this->acceptor_.close ();
      }
      return 0;
    }

    void add_node(CloudBusNodeController* c, GadgetronNodeInfo n)
    {
      mtx_.acquire();
      auto t = std::chrono::system_clock::from_time_t(n.last_recon);
      std::chrono::duration<double> time_since_last_recon =
	std::chrono::system_clock::now() - t;
      GDEBUG("Adding node: %s, %s, %d, (active reconstructions: %d, last recon %f s)\n",
	     n.uuid.c_str(), n.address.c_str(), n.port, n.active_reconstructions, time_since_last_recon.count());
      node_map_[c] = n;
      mtx_.release();
    }

    void delete_node(CloudBusNodeController* c)
    {
      mtx_.acquire();
      std::map<CloudBusNodeController*,GadgetronNodeInfo>::iterator it = node_map_.find(c);
      if (it != node_map_.end()) {
	GadgetronNodeInfo n = it->second;
	GDEBUG("Deleting node: %s, %s, %d\n", n.uuid.c_str(), n.address.c_str(), n.port);
	node_map_.erase(c);
      }
      mtx_.release();
    }

    void get_node_list(std::vector<GadgetronNodeInfo>& nl, CloudBusNodeController* exclude = 0)
    {
      mtx_.acquire();
      std::map<CloudBusNodeController*, GadgetronNodeInfo>::iterator it = node_map_.begin();
      nl.clear();
      while (it != node_map_.end()) {
	if (it->first != exclude) {
	  nl.push_back(it->second);
	}
	it++;
      }
      mtx_.release();
    }

  protected:
    ACE_SOCK_Acceptor acceptor_;
    std::map<CloudBusNodeController*, GadgetronNodeInfo> node_map_;
    ACE_Thread_Mutex mtx_;
  };



  class CloudBusNodeController
    : public ACE_Svc_Handler<ACE_SOCK_STREAM, ACE_MT_SYNCH>
  {
  public:
    CloudBusNodeController()
      : notifier_ (0, this, ACE_Event_Handler::WRITE_MASK)
    {
    }

    virtual ~CloudBusNodeController()
    {
      this->acceptor_->delete_node(this);
    }

    virtual int open (void) {
      this->notifier_.reactor (this->reactor ());
      this->msg_queue ()->notification_strategy (&this->notifier_);

      ACE_TCHAR peer_name[MAXHOSTNAMELENGTH];
      ACE_INET_Addr peer_addr;
      if (peer().get_remote_addr (peer_addr) == 0 &&
          peer_addr.addr_to_string(peer_name, MAXHOSTNAMELENGTH) == 0) {
	GDEBUG("Connection from %s\n", peer_name);
      }


      return this->reactor ()->register_handler(this,
						ACE_Event_Handler::READ_MASK);// | ACE_Event_Handler::WRITE_MASK);

    }

    virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE)
    {
      uint32_t msg_size;
      size_t recv_cnt;
      if ((recv_cnt = peer().recv_n (&msg_size, sizeof(uint32_t))) <= 0) {
	GERROR("Unable to read message message size\n");
	return -1;
      }

      char* buffer = new char[msg_size];

      if ((recv_cnt = peer().recv_n (buffer, msg_size)) <= 0) {
	GERROR("Unable to read message\n");
	delete [] buffer;
	return -1;
      }

      uint32_t msg_id = *((uint32_t*)buffer);
      GadgetronNodeInfo n;
      ACE_TCHAR peer_name[MAXHOSTNAMELENGTH];
      ACE_INET_Addr peer_addr;
      std::vector<GadgetronNodeInfo> nl;

      switch (msg_id) {
      case (GADGETRON_CLOUDBUS_NODE_INFO):

	deserialize(n,buffer+4,msg_size-4);

	//We will change the host name to the actually connected peer address since the relay host may not be able to resolve the host name
	if (peer().get_remote_addr (peer_addr) == 0) {
	  n.address = std::string(peer_addr.get_host_addr());
	}

	this->acceptor_->add_node(this,n);
	break;
      case (GADGETRON_CLOUDBUS_NODE_LIST_QUERY):
	{
	  //Get list of all nodes except myself
	  this->acceptor_->get_node_list(nl,this);
	  size_t buf_len = calculate_node_info_list_length(nl) + 4;
	  char* buffer = new char[buf_len+8];
	  *((uint32_t*)buffer) = buf_len+4;
	  *((uint32_t*)(buffer+4)) = GADGETRON_CLOUDBUS_NODE_LIST_REPLY;
	  try {
	    serialize(nl,buffer+8,buf_len);
	    this->peer().send_n(buffer,buf_len+8);
	  } catch (...) {
	    GERROR("Error serializing and sending node list\n");
	    throw;
	  }
	  break;
	}
      default:
	GERROR("Unknow message ID = %d\n", msg_id);
      }

      delete [] buffer;

      return 0;
    }

    virtual int handle_close (ACE_HANDLE handle,
			      ACE_Reactor_Mask mask)
    {
      GDEBUG("CloudBus connection connection closed\n");

      if (mask == ACE_Event_Handler::WRITE_MASK)
	return 0;

      this->stream_.close();

      mask = ACE_Event_Handler::ALL_EVENTS_MASK | ACE_Event_Handler::DONT_CALL;

      this->reactor ()->remove_handler (this, mask);

      //We are done with this controller.
      delete this;

      return 0;
    }

    virtual int output_ready(ACE_Message_Block* mb)
    {
      return 0;
    }

    void set_acceptor(CloudBusRelayAcceptor* a)
    {
      acceptor_ = a;
    }

  private:
    ACE_Reactor_Notification_Strategy notifier_;
    ACE_Stream<ACE_MT_SYNCH> stream_;
    CloudBusRelayAcceptor* acceptor_;
  };

  int CloudBusRelayAcceptor::handle_input (ACE_HANDLE fd)
  {
    CloudBusNodeController *controller;

    ACE_NEW_RETURN (controller, CloudBusNodeController, -1);

    controller->set_acceptor(this);

    if (this->acceptor_.accept (controller->peer ()) == -1) {
      GERROR("Failed to accept controller connection\n");
      delete controller;
      return -1;
    }

    controller->reactor (this->reactor ());
    if (controller->open () == -1)
      controller->handle_close (ACE_INVALID_HANDLE, 0);
    return 0;
  }


}


int main(int argc, char** argv)
{
  int port_no = 8002;
  unsigned short rest_port = 18002;
  
  if (argc > 3) {
    GERROR("Invalid number of arguments\n");
    GERROR("Usage: %s <port no> <rest port no>\n", argv[0]);
  }

  if (argc > 1) {
    port_no = std::atoi(argv[1]);
  }

  if (argc > 2) {
    rest_port = std::atoi(argv[2]);
  }

  
  ACE_INET_Addr port_to_listen (port_no);

  Gadgetron::CloudBusRelayAcceptor acceptor;

  acceptor.reactor (ACE_Reactor::instance ());
  if (acceptor.open (port_to_listen) == -1)
    return 1;

  Gadgetron::ReST::port_ = rest_port;
  
  Gadgetron::ReST::instance()->server()
    .route_dynamic("/info")([]()
			    {
			      return "Gadgetron CloudBusRelay\n";
			    });

  Gadgetron::ReST::instance()->server()
    .route_dynamic("/info/json")([&acceptor]()
				 {
				   std::vector<Gadgetron::GadgetronNodeInfo> nl;
				   acceptor.get_node_list(nl);
				   crow::json::wvalue out, nodes;
				   out["number_of_nodes"] = nl.size();
				   size_t idx = 0;
				   for (auto n: nl)
				   {
				     crow::json::wvalue node;
				     out["nodes"][idx]["uuid"] = n.uuid;
				     out["nodes"][idx]["address"] = n.address;
				     out["nodes"][idx]["port"] = n.port;
				     out["nodes"][idx]["rest_port"] = n.rest_port;
				     out["nodes"][idx]["compute_capability"] = n.compute_capability;
				     out["nodes"][idx]["active_reconstructions"] = n.active_reconstructions;
				     auto t = std::chrono::system_clock::from_time_t(n.last_recon);
				     std::chrono::duration<double> time_since_last_recon =
				       std::chrono::system_clock::now() - t;
				     out["nodes"][idx]["last_recon"] = time_since_last_recon.count();
				     idx++;
				   }
				   return out;
				 });

  ACE_Reactor::instance()->run_reactor_event_loop ();

  return 0;
}
