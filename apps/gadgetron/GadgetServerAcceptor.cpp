#include "GadgetServerAcceptor.h"
#include "GadgetStreamController.h"

using namespace Gadgetron;

GadgetServerAcceptor::~GadgetServerAcceptor ()
{
  this->handle_close (ACE_INVALID_HANDLE, 0);
}

int GadgetServerAcceptor::open (const ACE_INET_Addr &listen_addr)
{
  if (this->acceptor_.open (listen_addr, 1) == -1) {
    GERROR("error opening acceptor\n");
    return -1;
  }

  return this->reactor ()->register_handler(this, ACE_Event_Handler::ACCEPT_MASK);
}

int GadgetServerAcceptor::handle_input (ACE_HANDLE)
{
  GadgetStreamController *controller;

  ACE_NEW_RETURN (controller, GadgetStreamController, -1);

  auto_ptr<GadgetStreamController> p (controller);

  controller->set_global_gadget_parameters(global_gadget_parameters_);

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

int GadgetServerAcceptor::handle_close (ACE_HANDLE, ACE_Reactor_Mask)
{
  GDEBUG("GadgetServerAcceptor::handle_close\n");
  GDEBUG("Close Data Acceptor\n");

  if (this->acceptor_.get_handle () != ACE_INVALID_HANDLE) {
    ACE_Reactor_Mask m = 
      ACE_Event_Handler::ACCEPT_MASK | ACE_Event_Handler::DONT_CALL;
    this->reactor ()->remove_handler (this, m);
    this->acceptor_.close ();
  }
  return 0;
}
