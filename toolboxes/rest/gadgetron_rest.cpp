#include "gadgetron_rest.h"
#include <thread>

namespace Gadgetron
{
  namespace REST
  {
    Server<HTTP>* Server<HTTP>::instance_ = nullptr;
    unsigned int Server<HTTP>::port_ = 8080;
    size_t Server<HTTP>::threads_ = 1;
    
    Server<HTTP>* Server<HTTP>::instance()
    {
      if (!instance_) {
	instance_ = new Server<HTTP>(port_,threads_);
	Server<HTTP>* tmp = instance_;
	std::thread server_thread([tmp](){
	    tmp->start();
	});
	    
      }
      return instance_;
    }
  }
}
