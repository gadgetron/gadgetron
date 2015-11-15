#include "gadgetron_rest.h"

namespace Gadgetron
{
  ReST* Gadgetron::ReST::instance_ = nullptr;
  unsigned int Gadgetron::ReST::port_ = 9080;
  
  ReST* ReST::instance()
  {
    if (!instance_) {
      instance_ = new ReST();
	CROW_ROUTE(instance_->app_, "/")([]() {
	    return "<html><body><h1>GADGETRON</h1></body></html>\n";
	  });

	instance_->open();
      }
      return instance_;
    }

    void ReST::open()
    {
      Gadgetron::ReST* tmp = this;
      server_thread_ = std::thread([tmp](){
	  tmp->app_.port(port_)
	  .multithreaded()
	  .run();
	});
    }
}
