#pragma once

#include "Handlers.h"
#include "io/readers.h"
#include "Writers.h"
#include <tuple>
namespace Gadgetron::Server::Connection {
    template<class CONNECTION_TYPE>
    struct Connection {
        template<class ...ARGS>
        static inline auto process(ARGS&&... xs ) {

            auto connection = CONNECTION_TYPE(xs...);

            return connection.process();
        }
    };

    template<class HANDLER_FACTORY>
    inline void handle_input(std::istream& stream, HANDLER_FACTORY factory){
        bool closed = false;
        auto handlers = factory(closed);
        handlers[Handlers::CLOSE]    = std::make_unique<Handlers::CloseHandler>(closed);

        while (!closed) {
            auto id = Core::IO::read<uint16_t>(stream);

            GDEBUG_STREAM("Processing message with id: " << id);

            handlers.at(id)->handle(stream);
        }
    }

    inline void handle_output(std::shared_ptr<Core::InputChannel<Core::Message>> channel, std::ostream& stream, std::vector<std::unique_ptr<Core::Writer>> writers = {} ){
        using namespace Writers;
        writers.push_back(std::make_unique<TextWriter>());
        writers.push_back(std::make_unique<ResponseWriter>());

        for (auto message : *channel) {

            auto writer = std::find_if(writers.begin(), writers.end(),
                   [&](auto &writer) { return writer->accepts(*message); }
            );

            (*writer)->write(stream, std::move(message));
        }
    }


    template<class RETURN_TYPE>
    class BasicConnection {
    public:
        BasicConnection(std::iostream &stream, Gadgetron::Core::Context::Paths paths) :stream(stream),paths(paths){

         output_thread = std::thread(
                 [&](){handle_output(channel,stream);}
         );



        }
        RETURN_TYPE process();

        ~BasicConnection(){
            channel->close();
            output_thread.join();
        };

    private:

        std::thread output_thread;
        std::shared_ptr<Core::MessageChannel> channel = std::make_shared<Core::MessageChannel>();
        std::iostream& stream;
        const Gadgetron::Core::Context::Paths paths;
    };

}