
#include "ProtoConnection.h"

#include <typeindex>
#include <iostream>
#include <sstream>
#include <string>
#include <boost/dll/shared_library.hpp>
#include <boost/dll.hpp>
#include <boost/range/algorithm/transform.hpp>

#include "log.h"
#include "gadgetron_config.h"

#include "readers/Primitives.h"
#include "writers/ResponseWriter.h"
#include "Response.h"

#include "Reader.h"
#include "Writer.h"
#include "Server.h"
#include "Config.h"

#include "Builders.h"
#include "Handlers.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Core::Readers;
    using namespace Gadgetron::Server::Connection;
    using namespace Gadgetron::Server::Connection::Handlers;

    std::string read_filename_from_stream(std::istream &stream) {
        char buffer[1024];
        read_into(stream, buffer);
        return std::string(buffer);
    }

    class ConfigReferenceHandler : public Handler {
    public:
        ConfigReferenceHandler(
                std::function<void(Config config)> &escalate_callback,
                const Context::Paths &paths
        ) : callback(escalate_callback), paths(paths) {}

        void handle(std::istream &stream) override {
            boost::filesystem::path filename = paths.gadgetron_home / GADGETRON_CONFIG_PATH / read_filename_from_stream(stream);

            GDEBUG_STREAM("Reading config file: " << filename << std::endl);

            std::ifstream config_stream(filename.string());
            callback(parse_config(config_stream));
        }

    private:
        std::function<void(Config config)> &callback;
        const Context::Paths &paths;
    };

    class ConfigStringHandler : public Handler {
    public:
        explicit ConfigStringHandler(std::function<void(Config config)> &escalate_callback)
                : callback(escalate_callback) {}

        void handle(std::istream &stream) override {
            std::stringstream config_stream(read_string_from_stream<uint32_t>(stream));
            callback(parse_config(config_stream));
        }

    private:
        std::function<void(Config config)> &callback;
    };

    class CloseHandler : public Handler {
    public:
        explicit CloseHandler(bool &closed) : closed(closed) {}

        void handle(std::istream &stream) override {
            closed = true;
        }

    private:
        bool &closed;
    };

    class ReaderHandler : public Handler {
    public:
        ReaderHandler(std::unique_ptr<Reader> &&reader, std::shared_ptr<MessageChannel> channel)
                : reader(std::move(reader)), channel(std::move(channel)) {}

        void handle(std::istream &stream) override {
            channel->push_message(reader->read(stream));
        }

        std::unique_ptr<Reader> reader;
        std::shared_ptr<MessageChannel> channel;
    };
};



namespace Gadgetron::Server::Connection {

    void ProtoConnection::start() {
        auto self = shared_from_this();

        auto input_thread = std::thread ([=]() {
            self->process_input();
        });

        output_thread = std::thread([=]() {
            self->process_output();
        });

        input_thread.detach();
    };


    void ProtoConnection::process_input() {

        GDEBUG_STREAM("Input thread running.");

        bool closed = false;

        std::map<uint16_t, std::unique_ptr<Handler>> handlers{};

        std::function<void(Config)> escalate_callback = [&](const Config &config) {
            closed = true;
            this->escalate(config);
        };

        handlers[FILENAME] = std::make_unique<ConfigReferenceHandler>(escalate_callback, paths);
        handlers[CONFIG]   = std::make_unique<ConfigStringHandler>(escalate_callback);
        handlers[HEADER]   = std::make_unique<ErrorProducingHandler>("Received ISMRMRD header before config file.");
        handlers[CLOSE]    = std::make_unique<CloseHandler>(closed);
        handlers[QUERY]    = std::make_unique<QueryHandler>(channel);

        while (!closed) {
            auto id = read_t<uint16_t>(*stream);
            handlers.at(id)->handle(*stream);
        }
    }

    void ProtoConnection::process_output() {
        GDEBUG_STREAM("Output thread running.");

        InputChannel<Message>& input = *channel;
        for (auto message : input){
            GDEBUG("Hi! Listen!\n");
        }

    }

    void ProtoConnection::escalate(Config config) {
        channel->close();
        output_thread.join();


    }


}