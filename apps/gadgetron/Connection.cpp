
#include <typeindex>
#include <iostream>
#include <sstream>
#include <string>

#include "log.h"
#include "gadgetron_config.h"

#include "Connection.h"

#include "Stream.h"
#include "Server.h"
#include "Config.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Server;

    enum class message_id : uint16_t {
        FILENAME = 1,
        CONFIG   = 2,
        HEADER   = 3,
        CLOSE    = 4,
        TEXT     = 5,
        QUERY    = 6
    };

    class unexpected_message_type : public std::runtime_error {
    public:
        unexpected_message_type(uint16_t mid) : std::runtime_error(message(mid)) {}

    private:
        static std::string message(uint16_t mid) {
            std::stringstream message;
            message << "Unexpected message id: " << mid;
            return message.str();
        }
    };

    template <class T>
    void read_into(std::istream &stream, T &t) {
        stream.read(reinterpret_cast<char*>(&t), sizeof(t));
    }

    template <class T>
    T read_t(std::istream &stream) {
        T t;
        read_into(stream, t);
        return t;
    }

    std::string read_filename_from_stream(std::istream &stream) {
        char buffer[1024];
        read_into(stream, buffer);
        return std::string(buffer);
    }

    std::string read_string_from_stream(std::istream &stream) {

        uint32_t n = read_t<uint32_t>(stream);
        auto buffer = std::make_unique<char[]>(n);

        stream.read(buffer.get(), n);

        return std::string(buffer.get());
    }

    class ConfigFileHandler {
    public:
        ConfigFileHandler(std::promise<Gadgetron::Server::Config> &config_promise, const Context::Paths &paths)
                : promise_(config_promise), paths_(paths) {}

        void operator()(std::iostream &stream) {
            boost::filesystem::path filename = paths_.gadgetron_home / GADGETRON_CONFIG_PATH / read_filename_from_stream(stream);
            std::ifstream config_stream(filename.string());

            GDEBUG_STREAM("Reading config file: " << filename << std::endl);

            promise_.set_value(parse_config(config_stream));
        }

    private:
        std::promise<Config> &promise_;
        const Context::Paths &paths_;
    };

    class ConfigHandler {
    public:
        ConfigHandler(std::promise<Config> &config_promise) : promise_(config_promise){}

        void operator()(std::iostream &stream) {
            std::stringstream config_stream(read_string_from_stream(stream));

            promise_.set_value(parse_config(config_stream));
        }

    private:
        std::promise<Config> &promise_;
    };

    class HeaderHandler {
    public:
        HeaderHandler(std::promise<Context::Header> &header_promise) : promise_(header_promise) {}

        void operator()(std::iostream &stream) {
            std::string raw_header(read_string_from_stream(stream));

            ISMRMRD::IsmrmrdHeader header;
            ISMRMRD::deserialize(raw_header.c_str(), header);

            promise_.set_value(header);
        }

    private:
        std::promise<Context::Header> &promise_;
    };

    class CloseHandler {
    public:
        void operator()(std::iostream &stream) {
            throw std::runtime_error("Someone closed a thing.");
        }
    };

    class QueryHandler {
    public:
        void operator()(std::iostream &stream) {

        }
    };

    struct Chain {
        // Gadgetron::Core::Stream stream;
        std::map<message_id, std::shared_ptr<Reader>> readers;
    };

    Chain build_processing_chain(std::future<Config> config,
            std::future<Context::Header> header,
            const Context::Paths paths) {

        GINFO("Hello, I'm building a processing chain.\n");

        config.get();

        GINFO("Nice - a config! Thank you!\n");

        header.get();

        GINFO("Nice - a header! I can do things now!\n");

        return Chain();
    }

    void add_reader_to_handlers(std::map<message_id, std::function<void(std::iostream&)>> &map, message_id id, std::future<Chain> &chain) {

        GDEBUG_STREAM("Adding reader for message id: " << static_cast<int>(id) << std::endl);

        auto reader = chain.get().readers.at(id);

        map[id] = [=](std::iostream &stream) {
            auto message = reader->read(stream);
        };
    }
}



Connection::Connection(Context::Paths &paths, tcp::socket &socket)
    : stream_(std::move(socket)), paths_(paths) {}

void Connection::start() {

    auto self = shared_from_this();

    std::thread input_thread([=]() {
        self->process_input();
    });

    std::thread output_thread([=]() {
        self->process_output();
    });

    input_thread.detach();
    output_thread.detach();
}

void Connection::process_input() {

    GDEBUG_STREAM("Input thread running.");

    std::promise<Config> config_promise;
    std::promise<Context::Header> header_promise;

    std::future<Chain> chain = std::async(
            build_processing_chain,
            config_promise.get_future(),
            header_promise.get_future(),
            paths_
    );

    std::map<message_id, std::function<void(std::iostream&)>> handlers = {
        {message_id::FILENAME, ConfigFileHandler(config_promise, paths_)},
        {message_id::CONFIG, ConfigHandler(config_promise)},
        {message_id::HEADER, HeaderHandler(header_promise)},
        {message_id::CLOSE, CloseHandler()},
        {message_id::QUERY, QueryHandler()}
    };

    while (true) {
        message_id id = read_t<message_id>(stream_);

        GDEBUG_STREAM("Handling message with id: " << static_cast<int>(id) << std::endl);

        if (!handlers.count(id)) {
            add_reader_to_handlers(handlers, id, chain);
        }

        auto handler = handlers.at(id);
        handler(stream_);
    }
}

void Connection::process_output() {

    GDEBUG_STREAM("Output thread running.");


}

std::shared_ptr<Connection> Connection::create(Gadgetron::Core::Context::Paths &paths, tcp::socket &socket) {

    auto connection = std::shared_ptr<Connection>(new Connection(paths, socket));
    connection->start();

    return connection;
}
