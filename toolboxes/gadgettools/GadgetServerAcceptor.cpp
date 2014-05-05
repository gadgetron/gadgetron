#include "GadgetServerAcceptor.h"
#include "GadgetStreamController.h"

using namespace Gadgetron;

GadgetServerAcceptor::~GadgetServerAcceptor ()
{
  this->handle_close (ACE_INVALID_HANDLE, 0);
}

int GadgetServerAcceptor::open (const ACE_INET_Addr &listen_addr)
{
  if (this->acceptor_.open (listen_addr, 1) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
                       ACE_TEXT ("%p\n"),
                       ACE_TEXT ("acceptor.open")),
                      -1);
  return this->reactor ()->register_handler
    (this, ACE_Event_Handler::ACCEPT_MASK);
}

int GadgetServerAcceptor::handle_input (ACE_HANDLE)
{
  GadgetStreamController *controller;
  ACE_NEW_RETURN (controller, GadgetStreamController, -1);
  auto_ptr<GadgetStreamController> p (controller);

  if ( working_directory_.empty() )
  {
    #ifdef _WIN32
        working_directory_ = "c:\\temp\\gadgetron\\";
    #else
        working_directory_ = "/tmp/gadgetron/";
    #endif // _WIN32
  }

  controller->set_working_directory(working_directory_);

  if (this->acceptor_.accept (controller->peer ()) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
                       ACE_TEXT ("(%P|%t) %p\n"),
                       ACE_TEXT ("Failed to accept ")
                       ACE_TEXT ("controller connection")),
                      -1);
  p.release ();
  controller->reactor (this->reactor ());
  if (controller->open () == -1)
    controller->handle_close (ACE_INVALID_HANDLE, 0);
  return 0;
}

int GadgetServerAcceptor::handle_close (ACE_HANDLE, ACE_Reactor_Mask)
{
  ACE_DEBUG( (LM_DEBUG, 
	      ACE_TEXT("GadgetServerAcceptor::handle_close")) );
  
  GADGET_DEBUG1("Close Data Acceptor\n");

  if (this->acceptor_.get_handle () != ACE_INVALID_HANDLE) {
    ACE_Reactor_Mask m = 
      ACE_Event_Handler::ACCEPT_MASK | ACE_Event_Handler::DONT_CALL;
    this->reactor ()->remove_handler (this, m);
    this->acceptor_.close ();
  }
  return 0;
}
