#ifndef GADGETRON_REST_H
#define	GADGETRON_REST_H

#include <thread>

#include "server_http.hpp" //from https://github.com/eidheim/Simple-Web-Server

#include "gadgetron_rest_exports.h"

namespace Gadgetron {
    
  using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

  class EXPORTREST ReST
  {
  public:
    static ReST* instance();

    HttpServer* server()
    {
      return app_;
    }
    
    void restart();

    static unsigned int port_;
  protected:
    ReST() {};

    HttpServer* app_;

    std::thread server_thread_;
    static ReST* instance_;
  };
}
#endif	//GADGETRON_REST_H
