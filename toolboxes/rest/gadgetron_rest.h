#ifndef GADGETRON_REST_H
#define	GADGETRON_REST_H
#include <thread>
#include "crow_all.h"
#include "gadgetron_rest_exports.h"

namespace Gadgetron {
  class EXPORTREST ReST
  {
  public:
    static ReST* instance();
    crow::SimpleApp& server()
    {
      return app_;
    }
    
    static unsigned int port_;
  protected:
    ReST() {};
    void open();
    crow::SimpleApp app_;
    std::thread server_thread_;
    static ReST* instance_;
  };
}
#endif	//GADGETRON_REST_H
