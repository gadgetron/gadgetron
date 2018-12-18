
#include <memory>
#include <future>

#include "StreamConnection.h"

#include "Handlers.h"

#include "readers/Primitives.h"
#include "Channel.h"
#include "Context.h"


namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Readers;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Handlers;

    using Header = Gadgetron::Core::Context::Header;


    class HeaderHandler : public Handler {
    public:
        explicit HeaderHandler(std::function<void(Header)> callback) : callback{callback} {}

        void handle(std::istream &stream) override {
            std::string raw_header(read_string_from_stream<uint32_t>(stream));

            ISMRMRD::IsmrmrdHeader header;
            ISMRMRD::deserialize(raw_header.c_str(), header);
            callback(header);
        }

    private:
        std::function<void(Header)> callback;
    };


}

StreamConnection::StreamConnection(Core::Context context, Config config, std::unique_ptr<std::iostream> stream) {

}

void StreamConnection::start() {

}

void StreamConnection::process_input() {


    std::unordered_map<uint16_t, std::unique_ptr<Handler>> handlers;
    bool closed = false;

     while (!closed) {
        auto id = read_t<uint16_t>(*stream);
        handlers.at(id)->handle(*stream);
    }



}
