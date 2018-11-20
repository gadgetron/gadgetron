
#include <typeindex>
#include <iostream>
#include <sstream>
#include "Connection.h"

#include "log.h"
#include "gadgetron_config.h"

#include "Stream.h"
#include "Server.h"

namespace {

    using namespace Gadgetron::Core;

    enum class message_id : uint16_t {
        FILENAME = 1,
        CONFIG   = 2,
        HEADER   = 3,
        CLOSE    = 4,
        TEXT     = 5
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

    std::string read_config_string(std::istream &stream, Context::Paths &paths) {

        auto mid = read_t<message_id>(stream);

        if (mid == message_id::FILENAME) {
        }

        if (mid == message_id::CONFIG) {
            return read_string_from_stream(stream);
        }

        throw unexpected_message_type(static_cast<uint16_t>(mid));
    }

    std::string read_params_string(std::istream &stream) {

        auto mid = read_t<message_id>(stream);

        if (mid == message_id::HEADER) {
            return read_string_from_stream(stream);
        }

        throw unexpected_message_type(static_cast<uint16_t>(mid));
    }
}

class ConfigFileHandler {
public:
    ConfigFileHandler(std::promise<Context::Config> &config_promise, const Context::Paths &paths) : paths_(paths) {}

    void operator()(std::iostream &stream) {
        boost::filesystem::path filename = paths_.gadgetron_home / GADGETRON_CONFIG_PATH / read_filename_from_stream(stream);

        std::ifstream config_stream(filename.string());

        std::string config((std::istream_iterator<char>(config_stream)),
                            std::istream_iterator<char>());

        GINFO_STREAM(config);
    }

private:
    const Context::Paths &paths_;
};

class ConfigHandler {
public:
    ConfigHandler(std::promise<Context::Config> &config_promise) {}

    void operator()(std::iostream &stream) {

    }
};

class HeaderHandler {
public:
    HeaderHandler(std::promise<Context::Header> &header_promise) {}

    void operator()(std::iostream &stream) {

    }
};

class CloseHandler {
public:
    CloseHandler(std::promise<Context::Header> &header_promise) {}

    void operator()(std::iostream &stream) {

    }
};

Connection::Connection(Gadgetron::Core::Context::Paths &paths, tcp::socket &socket)
    : stream_(std::move(socket)), paths_(paths) {}

void Connection::start() {

    auto self = shared_from_this();

    std::thread connection_thread([=]() {
        self->process_input();
    });

    connection_thread.detach();
}

void Connection::process_input() {

    GDEBUG_STREAM("Connection thread running.");

    Context::Config config = read_config(stream_);
    Context::Header header = read_header(stream_);

    struct {
        std::map<message_id, Reader> readers;
        std::map<std::type_index, Writer> writers;
        InputChannel<Message> input;
        OutputChannel output;
    } stuff = build_streams_and_stuff(config, header);

    while(true) {
        message_id mid;
        read_into(stream_, mid);

        Reader &reader = stuff.readers.at(mid);
        stuff.input.put(reader.read(stream_));
    }
}

std::shared_ptr<Connection> Connection::create(Gadgetron::Core::Context::Paths &paths, tcp::socket &socket) {

    auto connection = std::make_shared<Connection>(paths, socket);
    connection->start();

    return connection;
}
