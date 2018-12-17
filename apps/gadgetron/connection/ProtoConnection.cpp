
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

//void start_stream(const Config &config, std::future<Header> header_future) {
//
//    Context::Header header = header_future.get();
//    Context context{header, paths};
//
//    auto stream = builder.build_stream(config.stream, context);
//
//    stream->process(channels.input, channels.output);
//}
//
//void initialize_readers(const Config &config) {
//
//    std::map<uint16_t, std::unique_ptr<Handler>> handlers;
//
//    auto readers = builder.build_readers(config.readers);
//
//    for (auto &reader_pair : readers) {
//        handlers.emplace(reader_pair.first,
//                         std::make_unique<ReaderHandler>(std::move(reader_pair.second), channels.input));
//    }
//    this->promises.readers.set_value(std::move(handlers));
//}
//
//void initialize_writers(const Config &config) {
//    this->promises.writers.set_value(builder.build_writers(config.writers));
//}

namespace Gadgetron::Server::Connection {

    void ProtoConnection::start() {
        auto self = shared_from_this();

        std::thread input_thread([=]() {
            self->process_input();
        });

        std::thread output_thread([=]() {
            self->process_output();
        });

        input_thread.detach();
        output_thread.detach();
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

//        auto writer_future = this->promises.writers.get_future();
//        auto writers = writer_future.get();
//
////        auto writers = std::vector<std::unique_ptr<Writer>>();
////        writers.push_back(std::make_unique<Writers::ResponseWriter>());
//
//        std::shared_ptr<InputChannel<Message>> output = this->channels.output;
//
////        for (std::unique_ptr<Message> message : *output) {
//        try {
//            while (true) {
//                auto message = output->pop();
//                GDEBUG_STREAM("Ptr " << message.get() << std::endl);
//                GDEBUG_STREAM("Writer got a: " << typeid(*message).name() << std::endl);
//                auto writer = std::find_if(writers.begin(), writers.end(),
//                                           [&](auto &writer) { return writer->accepts(*message); }
//                );
//
//
//                if (writer != writers.end()) {
//                    (*writer)->write(*stream, std::move(message));
//                }
//            }
//        } catch (ChannelClosedError err) {
//
//        }
    }
}