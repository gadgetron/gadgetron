
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
        explicit HeaderHandler(std::promise<Header> &header_promise) : promise(header_promise) {}

        void handle(std::istream &stream) override {
            std::string raw_header(read_string_from_stream<uint32_t>(stream));

            ISMRMRD::IsmrmrdHeader header;
            ISMRMRD::deserialize(raw_header.c_str(), header);

            promise.set_value(header);
        }

    private:
        std::promise<Header> &promise;
    };

}