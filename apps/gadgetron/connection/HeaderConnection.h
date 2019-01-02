#pragma once

#include "Connection.h"
#include "Config.h"
#include "Channel.h"
#include "Context.h"


namespace Gadgetron::Server::Connection {

    class HeaderConnection : public Connection {
    public:
        using Context = Gadgetron::Core::Context;
        using Header = Context::Header;

        static void process(
                std::iostream &stream,
                const Context::Paths &paths,
                const Config &config,
                ErrorHandler &error_handler
        );

    protected:
        HeaderConnection(std::iostream &stream, Gadgetron::Core::Context::Paths paths);

        std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(std::function<void()> close) override;

        std::promise<boost::optional<Header>> promise;
        const Gadgetron::Core::Context::Paths paths;
    };
}
