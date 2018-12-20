#ifndef GADGETRON_CONNECTION_H
#define GADGETRON_CONNECTION_H

#include <boost/asio.hpp>

#include "connection/Handlers.h"

#include "Writer.h"
#include "Context.h"
#include "Channel.h"

namespace Gadgetron::Server::Connection {
    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream);


    class Connection {
    public:
        using MessageChannel = Gadgetron::Core::MessageChannel;
        using Handler = Gadgetron::Server::Connection::Handlers::Handler;
        using Writer = Gadgetron::Core::Writer;

    protected:
        virtual std::map<uint16_t, std::unique_ptr<Handler>> prepare_handlers(bool &closed) = 0;
        virtual std::vector<std::unique_ptr<Writer>> prepare_writers();

        void start();
        void join();

        explicit Connection(std::iostream &stream);
        virtual ~Connection() = default;

        void process_input();
        void process_output();

        std::iostream &stream;

        struct {
            std::shared_ptr<MessageChannel> input, output;
        } channels;

        struct {
            std::thread input, output;
        } threads;
    };
}


#endif //GADGETRON_CONNECTION_H
