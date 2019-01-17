
#pragma once

#include <memory>

#include "connection/Loader.h"
#include "connection/stream/Processable.h"

#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Server::Connection::Stream {

    class Stream : public Processable {
    public:
        const std::string key;

        Stream(const Config::Stream &, const Core::Context &, Loader &);

        void process(
                std::shared_ptr<Core::Channel> input,
                std::shared_ptr<Core::Channel> output,
                ErrorHandler &error_handler
        ) override;

    private:
        std::vector<std::unique_ptr<Processable>> nodes;
    };
}
