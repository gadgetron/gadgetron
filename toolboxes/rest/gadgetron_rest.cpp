#include "gadgetron_rest.h"

namespace Gadgetron
{
  namespace ReST
  {
    GadgetronReST* GadgetronReST::instance_ = nullptr;
    unsigned int GadgetronReST::port_ = 9080;
    size_t GadgetronReST::threads_ = 1;
    
    GadgetronReST* GadgetronReST::instance()
    {
      if (!instance_) {
	instance_ = new GadgetronReST(port_,threads_);

	instance_->default_resource["GET"]=[](GadgetronReST::Response& response, std::shared_ptr<GadgetronReST::Request> request)
	  {
	    std::string content = "<html><body><h1>GADGETRON</h1></body></html>\n";
	    response << "HTTP/1.1 200 OK\r\nContent-Length: " << content.length() << "\r\n\r\n" << content;
	  };
	
	instance_->open();
      }
      return instance_;
    }

    void GadgetronReST::open()
    {
      GadgetronReST* tmp = this;
      server_thread_ = std::thread([tmp](){
	  tmp->start();
	});
    }
  }
}
