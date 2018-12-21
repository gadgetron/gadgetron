
#include <readers/Primitives.h>
#include "ConfigConnection.h"

#include "Context.h"

#include "Handlers.h"
#include "Writers.h"
#include "Config.h"
#include "Connection_common.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Readers;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Writers;
    using namespace Gadgetron::Server::Connection::Handlers;

    using Header = Gadgetron::Core::Context::Header;

    class HeaderHandler : public Handler {
    public:
        explicit HeaderHandler(
                std::function<void(Header)> header_callback
        ) : header_callback(std::move(header_callback)) {}

        void handle(std::istream &stream) override {
            std::string raw_header(read_string_from_stream<uint32_t>(stream));

            ISMRMRD::IsmrmrdHeader header{};
            ISMRMRD::deserialize(raw_header.c_str(), header);

            header_callback(header);
        }

    private:
        std::function<void(Header)> header_callback;
    };
}

#define CONFIG_ERROR "Received second config file. Only one config allowed."

namespace Gadgetron::Server::Connection {

        using Header = Core::Context::Header ;
        template<> Header ConfigConnection::process() {

        Header header{};

        auto factory = [&](auto& closed ) {

                auto header_callback = [&](Header input_header) { header = input_header; closed = true;};

                std::map<uint16_t, std::unique_ptr<Handler>> handlers{};
                handlers[FILENAME] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
                handlers[CONFIG] = std::make_unique<ErrorProducingHandler>(CONFIG_ERROR);
                handlers[QUERY] = std::make_unique<QueryHandler>(channel);

                handlers[HEADER] = std::make_unique<HeaderHandler>(header_callback);
                return handlers;
            };

        handle_input(stream,factory);

        return header;
    }

    template class BasicConnection<Header>;



}