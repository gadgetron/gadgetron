#ifndef GADGETRON_SERVER_H
#define GADGETRON_SERVER_H


#include <boost/filesystem/path.hpp>
#include "storage_server.h"

namespace Gadgetron::Server {

  struct Settings {
    Core::Context::Paths paths;
    std::uint16_t port;
    StorageServer::Address storage_address;
  };

    class Server {
    public:
        explicit Server(const Settings& settings);
        [[noreturn]] void serve();

    private:
      Settings settings;

    };
}

#endif //GADGETRON_SERVER_H
