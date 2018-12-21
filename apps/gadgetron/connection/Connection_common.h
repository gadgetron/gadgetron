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
        BasicConnection(std::iostream& stream, Gadgetron::Core::Context::Paths paths) :stream(stream),paths(paths){

        }

        RETURN_TYPE process(ErrorHandler& handler){

            auto channel = std::make_shared<MessageChannel>();

            auto output_thread = std::thread(
                    [&]() {handler.handle("Connection Output Thread",
                            [&](){
                                handle_output(channel, stream);
                            })
                    });

            return process_input();

        }

    private:

        RETURN_TYPE process_input();

        std::iostream& stream;
        const Gadgetron::Core::Context::Paths paths;
    };

}