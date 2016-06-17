#ifndef GADGETRON_REST_H
#define	GADGETRON_REST_H

#include <thread>

#include "crow/crow.h"
#include "crow/json.h"

// require these undef to avoid visual studio conflicts for ERROR and DELETE
// TODO: remove these once updated to visual studio 2015
#ifdef WIN32
    #ifdef ERROR
        #undef ERROR
    #endif // ERROR

    #ifdef DELETE
        #undef DELETE
    #endif // DELETE
#endif // WIN32

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
