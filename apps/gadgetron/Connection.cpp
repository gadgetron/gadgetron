
#include <memory>
#include <connection/StreamConnection.h>

#include "connection/ProtoConnection.h"
#include "connection/ConfigConnection.h"
#include "connection/Writers.h"

#include "Connection.h"
#include "connection/Connection_common.h"


#include "readers/Primitives.h"
#include "log.h"

using namespace boost::asio;

using namespace Gadgetron::Core;
using namespace Gadgetron::Core::Readers;

using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Writers;

namespace {

    void send_errors(std::iostream &stream) {}

    void send_close(std::iostream &stream) {
        uint16_t close = 4;
        stream.write(reinterpret_cast<char *>(&close), sizeof(close));
    }

}

void Gadgetron::Server::Connection::handle(const Gadgetron::Core::Context::Paths &paths,
                                           std::unique_ptr<std::iostream> stream) {

    auto config = process<ProtoConnection, decltype(*stream), decltype(paths)>(*stream,paths);

    if (config){

        auto header = process<ConfigConnection, decltype(*stream), decltype(paths)>(*stream,paths);

        auto context = Core::Context{header,paths};
        process<StreamConnection,decltype(*stream),decltype(context),decltype(*config)>(*stream,context,*config);
        GDEBUG("WEEEE");
    }

    send_errors(*stream);
    send_close(*stream);
}
