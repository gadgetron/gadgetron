
#include <memory>

#include "connection/VoidConnection.h"
#include "connection/ConfigConnection.h"
#include "connection/HeaderConnection.h"
#include "connection/StreamConnection.h"
#include "connection/Writers.h"
#include "connection/Core.h"

#include "Connection.h"


#include "readers/Primitives.h"
#include "log.h"

using namespace boost::asio;

using namespace Gadgetron::Core;
using namespace Gadgetron::Core::Readers;

using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Writers;


namespace Gadgetron::Server::Connection {

    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {
        auto thread = std::thread(handle_connection, std::move(stream), paths);
        thread.detach();
    }
}

