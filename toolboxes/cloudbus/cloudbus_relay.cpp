#include "ace/Reactor.h"
#include "ace/SOCK_Acceptor.h"
#include "ace/Svc_Handler.h"
#include "ace/SOCK_Stream.h"
#include "ace/Reactor_Notification_Strategy.h"
#include "ace/Stream.h"

#include "log.h"

namespace Gadgetron{

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
    }
    
    virtual int open (void) {
      this->notifier_.reactor (this->reactor ());
      this->msg_queue ()->notification_strategy (&this->notifier_);

      ACE_TCHAR peer_name[MAXHOSTNAMELEN];
      ACE_INET_Addr peer_addr;
      if (peer().get_remote_addr (peer_addr) == 0 &&
	  peer_addr.addr_to_string (peer_name, MAXHOSTNAMELEN) == 0) {
	GINFO("Connection from %s\n", peer_name);
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
      
      GDEBUG("Received message of size: %d\n", msg_size);
      delete [] buffer;

      return 0;
    }

    virtual int handle_close (ACE_HANDLE handle,
			      ACE_Reactor_Mask mask)
    {
      GDEBUG("Cloud bus connection closed\n");
      
      if (mask == ACE_Event_Handler::WRITE_MASK)
	return 0;
      
      this->stream_.close();
      
      mask = ACE_Event_Handler::ALL_EVENTS_MASK |
	ACE_Event_Handler::DONT_CALL;
  
      this->reactor ()->remove_handler (this, mask);
    }
    
    virtual int output_ready(ACE_Message_Block* mb) 
    {
      return 0;
    }

  private:
    ACE_Reactor_Notification_Strategy notifier_;
    ACE_Stream<ACE_MT_SYNCH> stream_;
    
  };



  class CloudBusRelayAcceptor : public ACE_Event_Handler
  {
  public:
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
    
    virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE)
    {
      CloudBusNodeController *controller;

      ACE_NEW_RETURN (controller, CloudBusNodeController, -1);
      
      auto_ptr<CloudBusNodeController> p (controller);
      
      if (this->acceptor_.accept (controller->peer ()) == -1) {
	GERROR("Failed to accept controller connection\n"); 
	return -1;
      }
      
      p.release ();
      controller->reactor (this->reactor ());
      if (controller->open () == -1)
	controller->handle_close (ACE_INVALID_HANDLE, 0);
      return 0;

    }
    
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
    
  protected:
    ACE_SOCK_Acceptor acceptor_;
  };
}


int main(int argc, char** argv)
{
  GDEBUG("CloudBus Relay Staring\n");

  const int port_no = 8002;

  ACE_INET_Addr port_to_listen (port_no);

  Gadgetron::CloudBusRelayAcceptor acceptor;

  acceptor.reactor (ACE_Reactor::instance ());
  if (acceptor.open (port_to_listen) == -1)
    return 1;

  ACE_Reactor::instance()->run_reactor_event_loop ();

  return 0;
}
