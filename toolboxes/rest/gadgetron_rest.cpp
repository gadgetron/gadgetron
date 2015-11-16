#include <chrono>
#include "gadgetron_rest.h"

namespace Gadgetron
{
  ReST* Gadgetron::ReST::instance_ = nullptr;
  unsigned int Gadgetron::ReST::port_ = 9080;
  
  ReST* ReST::instance()
  {
    if (!instance_) {
      instance_ = new ReST();
      instance_->app_.route_dynamic("/")([]() {
	  return "<html><body><h1>GADGETRON</h1></body></html>\n";
	});

	instance_->open();
      }
      return instance_;
    }

    void ReST::open()
    {
      crow::logger::setLogLevel(crow::LogLevel::ERROR);
      Gadgetron::ReST* tmp = this;
      server_thread_ = std::thread([tmp](){
	  tmp->app_.port(port_).run();
	});

      //TODO: We should get rid of this.
      //Not sure what the deal is, but we get a crash here if we return the instance before it
      //is properly started. 
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}
