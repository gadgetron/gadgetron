#pragma once

#include "Channel.h"
#include "Context.h"

#include "connection/Core.h"

namespace Gadgetron::Server::Connection::Stream {

    class DecoratedErrorHandler : public ErrorHandler {
    public:
        ErrorHandler &handler;
        std::string decorator;

        DecoratedErrorHandler(ErrorHandler &handler, std::string decorator)
                : handler(handler), decorator(std::move(decorator)) {}

        void handle(const std::string &location, std::function<void()> function) override {
            handler.handle(decorator + "/" + location, function);
        }
    };



    class Processable {
    public:
        virtual ~Processable() = default;

        virtual void process(
                std::shared_ptr<Core::Channel> input,
                std::shared_ptr<Core::Channel> output,
                ErrorHandler &error_handler
        ) = 0;
    };
}