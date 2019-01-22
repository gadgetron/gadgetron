#pragma once

#include <thread>
#include <memory>
#include <functional>

#include "CloseGuard.h"

#include "io/primitives.h"
#include "Writer.h"
#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Server::Connection {

    class ErrorHandler {
    public:
        virtual ~ErrorHandler() = default;
        virtual void handle(const std::string &location, std::function<void()> function) = 0;

        template<class F, class... ARGS>
        std::thread run(const std::string &location, F fn, ARGS... args) {
            return std::thread(
                    [=](auto... args) {
                        auto decorated = [&]() { fn(args...); };
                        this->handle(location, decorated);
                    },
                    args...
            );
        }
    };

    template<class F>
    void process_input(std::iostream &stream, std::shared_ptr<Core::MessageChannel> input, F handler_factory) {

        CloseGuard closer{input};

        bool closed = false;
        auto handlers = handler_factory([&]() { closed = true; });

        while (!closed) {
            auto id = Core::IO::read<uint16_t>(stream);
            handlers.at(id)->handle(stream);
        }
    }

    template<class F>
    void process_output(std::iostream &stream, std::shared_ptr<Core::MessageChannel> output, F writer_factory) {

        CloseGuard closer{output};

        auto writers = writer_factory();

        Core::InputChannel &messages = *output;
        for (auto message : messages) {

            auto writer = std::find_if(writers.begin(), writers.end(),
                    [&](auto &writer) { return writer->accepts(message); }
            );

            if (writer != writers.end()) {
                (*writer)->write(stream, message);
            }
        }
    }

    std::vector<std::unique_ptr<Core::Writer>> default_writers();

    template<class F>
    std::thread start_input_thread(
            std::iostream &stream,
            std::shared_ptr<Core::MessageChannel> channel,
            F handler_factory,
            ErrorHandler &error_handler
    ) {
        return error_handler.run(
                "Connection Input Thread",
                [&](auto... args) { process_input(stream, args...); },
                channel, handler_factory
        );
    }

    template<class F>
    std::thread start_output_thread(
            std::iostream &stream,
            std::shared_ptr<Core::MessageChannel> channel,
            F writer_factory,
            ErrorHandler &error_handler
    ) {
        return error_handler.run(
                "Connection Output Thread",
                [&](auto... args) { process_output(stream, args...); },
                channel, writer_factory
        );
    }

    void handle_connection(std::unique_ptr<std::iostream> stream, Core::Context::Paths paths);
}