#pragma once
#include "Message.h"
#include "PureGadget.h"
#include "connection/Loader.h"
#include "connection/config/Config.h"

namespace Gadgetron::Server::Connection::Nodes {
    class PureStream {
    public:
        PureStream(const Config::PureStream&, const Core::Context&, Loader&);
        Core::Message process_function(Core::Message) const;

    private:
        const std::vector<std::unique_ptr<Core::GenericPureGadget>> pure_gadgets;
    };
}
