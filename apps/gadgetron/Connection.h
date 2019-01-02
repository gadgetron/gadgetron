#ifndef GADGETRON_CONNECTION_H
#define GADGETRON_CONNECTION_H

#include <boost/asio.hpp>

#include "connection/Handlers.h"

#include "Writer.h"
#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Connection {
    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream);

    class ErrorHandler {
    public:
        virtual ~ErrorHandler() = default;
        virtual void handle(const std::string &location, std::function<void()> function) = 0;

        template<class F, class... ARGS>
        std::thread run(const std::string &location, F fn, ARGS... args) {
            auto decorated = [=]() { fn(std::forward(args)...); };
            return std::thread(
                    [=]() {
                        this->handle(location, decorated);
                    }
            );
        }
    };

    class Connection {
    public:
        using MessageChannel = Gadgetron::Core::MessageChannel;
        using Writer = Gadgetron::Core::Writer;
        using Handler = Handlers::Handler;

    protected:
        explicit Connection(std::iostream &stream);
        virtual ~Connection() = default;

        virtual std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(std::function<void()> close) = 0;
        virtual std::vector<std::unique_ptr<Writer>> prepare_writers();

        void process_input();
        void process_output();

        std::iostream &stream;

        struct {
            std::shared_ptr<MessageChannel> input, output;
        } channels;
    };
}


#endif //GADGETRON_CONNECTION_H
