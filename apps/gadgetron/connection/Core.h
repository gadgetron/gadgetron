#pragma once

#include <thread>
#include <memory>
#include <functional>


#include "io/primitives.h"
#include "Writer.h"
#include "Channel.h"
#include "Context.h"
#include "Server.h"

namespace Gadgetron::Server::Connection {

    class ErrorReporter {
    public:
        virtual void operator()(const std::string&, const std::string&) = 0;
        virtual ~ErrorReporter() = default;
    };

    class ErrorHandler {
    public:

        ErrorHandler( ErrorReporter &callback, std::string location)
                : location{
                std::move(location)}, push_error{callback} {}

        ErrorHandler( const ErrorHandler &parent,const std::string &new_location) : push_error{parent.push_error}, location{
                add_location(parent.location, new_location)} {

        }

        ErrorHandler(const ErrorHandler &parent) = default;


#if defined(NDEBUG)

        template<class F, class...ARGS>
        void handle(F fn, ARGS &&... args) {
            // When debugging, it is useful to have all exceptions bubble up to the
            // debugger. To enable this, we sabotage the error handler on debug builds.
                try {
                    fn(std::forward<ARGS>(args)...);
                }
                catch (const Core::ChannelClosed &e) {
                    // Ignored.
                }
                catch (const std::exception &e) {
                    push_error(location, e.what());
                }
                catch (...) {
                    push_error(location, "Unknown error.");
                }
        }
#else

        template<class F, class...ARGS>
        void handle(F fn, ARGS &&... args) {
            try {
                fn(std::forward<ARGS>(args)...);
            }
            catch (const Core::ChannelClosed &e) {
                // Ignored.
            }
        }

#endif

        template<class F, class... ARGS>
        std::thread run(F fn, ARGS &&... args) {
            return std::thread(
                    []( auto handler, auto fn, auto &&... iargs) {
                        handler.handle(fn, std::forward<ARGS>(iargs)...);
                    },
                    *this,
                    std::forward<F>(fn),
                    std::forward<ARGS>(args)...
            );
        }

    private:

        static std::string add_location(const std::string &old_location, const std::string &new_location) {
            return old_location + "/" + new_location;
        }

        const std::string location;
        ErrorReporter &push_error;
    };

    template<class F>
    void process_input(std::iostream &stream, Core::OutputChannel channel, F handler_factory) {

        bool closed = false;
        auto handlers = handler_factory([&]() { closed = true; });

        while (!closed) {
            auto id = Core::IO::read<uint16_t>(stream);
            handlers.at(id)->handle(stream, channel);
        }
    }

    template<class F>
    void process_output(std::iostream &stream, Core::GenericInputChannel messages, F writer_factory) {

        auto writers = writer_factory();

        for (auto message : messages) {

            auto writer = std::find_if(writers.begin(), writers.end(),
                                       [&](auto &writer) { return writer->accepts(message); }
            );

            if (writer != writers.end()) {
                (*writer)->write(stream, std::move(message));
            }
        }
    }

    std::vector<std::unique_ptr<Core::Writer>> default_writers();

    template<class F>
    std::thread start_input_thread(
            std::iostream &stream,
            Core::OutputChannel channel,
            F handler_factory,
            ErrorHandler &error_handler
    ) {

        return ErrorHandler(error_handler,"Connection Input Thread").run(
                [&stream](auto c, auto h) { process_input(stream, std::move(c), h); },
                std::move(channel), handler_factory);
    }

    template<class F>
    std::thread start_output_thread(
            std::iostream &stream,
            Core::GenericInputChannel channel,
            F writer_factory,
            ErrorHandler &error_handler
    ) {
        return ErrorHandler(error_handler,"Connection Output Thread").run(
                [&stream](auto c, auto w) { process_output(stream, std::move(c), w); },
                std::move(channel), writer_factory
        );
    }

    void handle_connection(std::unique_ptr<std::iostream> stream, Settings settings);
}