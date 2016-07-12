#include "gadgetron_rest.h"
#include "GadgetServerAcceptor.h"
#include "GadgetStreamController.h"
#include "CloudBus.h"

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

  {
      std::lock_guard<std::mutex> guard(acceptor_mtx_);
      is_listening_ = true;
  }
  
  //Register a way to close the Acceptor via the ReST API
  Gadgetron::ReST::instance()->server()
      ->resource["/acceptor/close"]["GET"]=[this](std::shared_ptr<Gadgetron::HttpServer::Response> response, 
                                    std::shared_ptr<Gadgetron::HttpServer::Request> request) 
      {
          this->close();
          auto it = this->global_gadget_parameters_.find("using_cloudbus");
          if (it != this->global_gadget_parameters_.end() && it->second == std::string("true")) {
              CloudBus::set_relay_port(0);
              CloudBus::instance()->close();
          }

          std::stringstream content_stream;
          content_stream << "Acceptor closed\n";
          
          //find length of content_stream (length received using content_stream.tellp()) 
          content_stream.seekp(0, std::ios::end);
          
          *response <<  "HTTP/1.1 200 OK\r\nContent-Length: " << content_stream.tellp() << "\r\n\r\n" << content_stream.rdbuf();
      };


  //Register a way to get the port number if it is listening.
  Gadgetron::ReST::instance()->server()
      ->resource["/info/port"]["GET"]=[this,listen_addr](std::shared_ptr<Gadgetron::HttpServer::Response> response, 
                                         std::shared_ptr<Gadgetron::HttpServer::Request> request) 
      {

          std::stringstream content_stream;

          if (this->is_listening_) {
              std::stringstream ss;
              content_stream << listen_addr.get_port_number();
              //find length of content_stream (length received using content_stream.tellp()) 
              content_stream.seekp(0, std::ios::end);   
              *response <<  "HTTP/1.1 200 OK\r\nContent-Length: " << content_stream.tellp() << "\r\n\r\n" << content_stream.rdbuf();
              return;
          }
          
          content_stream << "Port not available";
          
          //find length of content_stream (length received using content_stream.tellp()) 
          content_stream.seekp(0, std::ios::end);
          
          *response <<  "HTTP/1.1 500 Internal Server Error\r\nContent-Length: " << content_stream.tellp() << "\r\n\r\n" << content_stream.rdbuf();
      };
  

  //Make sure latest resources load
  Gadgetron::ReST::instance()->restart();

  return this->reactor ()->register_handler(this, ACE_Event_Handler::ACCEPT_MASK);
}

int GadgetServerAcceptor::handle_input (ACE_HANDLE)
{
  GadgetStreamController *controller;

  ACE_NEW_RETURN (controller, GadgetStreamController, -1);


  controller->set_global_gadget_parameters(global_gadget_parameters_);

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
