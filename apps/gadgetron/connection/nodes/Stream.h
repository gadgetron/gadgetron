
#pragma once

#include <memory>

#include "connection/Loader.h"
#include "connection/config/Config.h"

#include "connection/core/Processable.h"

#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Server::Connection {
    class Loader;
}

namespace Gadgetron::Server::Connection::Nodes {

    class Stream : public Processable {
    public:
        const std::string key;
        Stream(const Config::Stream &, const Core::StreamContext &, Loader &);

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
