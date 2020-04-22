
#pragma once

#include <memory>

#include "connection/Config.h"
#include "connection/Loader.h"
#include "connection/stream/Processable.h"

#include "Channel.h"
#include "Context.h"
#include "Connection.h"

namespace Gadgetron::Server::Connection {
    class Loader;
}

namespace Gadgetron::Server::Connection::Stream {

    class Stream : public Processable {
    public:
        const std::string key;
        Stream(const Config::Stream &, const StreamContext &, Loader &);

        void process(
                Core::GenericInputChannel input,
                Core::OutputChannel output,
                ErrorHandler &
        ) override;

        bool empty() const;
        const std::string &name() override;

    private:
        std::vector<std::shared_ptr<Processable>> nodes;
    };
}

