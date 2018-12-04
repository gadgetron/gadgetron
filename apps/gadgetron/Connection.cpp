
#include <typeindex>
#include <iostream>
#include <sstream>
#include <string>
#include <boost/dll/shared_library.hpp>
#include <boost/dll.hpp>

#include "log.h"
#include "gadgetron_config.h"

#include "Connection.h"

#include "Builders.h"
#include "Stream.h"
#include "Server.h"
#include "Config.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Server;

    using Header = Gadgetron::Core::Context::Header;

    enum message_id : uint16_t {
        FILENAME    = 1,
        CONFIG      = 2,
        HEADER      = 3,
        CLOSE       = 4,
        TEXT        = 5,
        QUERY       = 6
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

    class Handler {
    public:
        virtual void handle(std::iostream &stream) = 0;
    };

    class ConfigFileHandler : public Handler {
    public:
        ConfigFileHandler(std::promise<Gadgetron::Server::Config> &config_promise, const Context::Paths &paths)
                : promise(config_promise), paths(paths) {}

        void handle(std::iostream &stream) override {
            boost::filesystem::path filename = paths.gadgetron_home / GADGETRON_CONFIG_PATH / read_filename_from_stream(stream);
            std::ifstream config_stream(filename.string());

            GDEBUG_STREAM("Reading config file: " << filename << std::endl);

            promise.set_value(parse_config(config_stream));
        }

    private:
        std::promise<Config> &promise;
        const Context::Paths &paths;
    };

    class ConfigHandler : public Handler {
    public:
        ConfigHandler(std::promise<Config> &config_promise) : promise(config_promise){}

        void handle(std::iostream &stream) override {
            std::stringstream config_stream(read_string_from_stream(stream));

            promise.set_value(parse_config(config_stream));
        }

    private:
        std::promise<Config> &promise;
    };

    class HeaderHandler : public Handler {
    public:
        HeaderHandler(std::promise<Header> &header_promise) : promise(header_promise) {}

        void handle(std::iostream &stream) override {
            std::string raw_header(read_string_from_stream(stream));

            ISMRMRD::IsmrmrdHeader header;
            ISMRMRD::deserialize(raw_header.c_str(), header);

            promise.set_value(header);
        }

    private:
        std::promise<Header> &promise;
    };

    class CloseHandler : public Handler {
    public:
        CloseHandler(bool &closed) : closed(closed) {}

        void handle(std::iostream &stream) override {
            closed = true;
        }
    private:
        bool &closed;
    };

    class QueryHandler : public Handler {
    public:
        void handle(std::iostream &stream) override {

        }
    };

    class ReaderHandler : public Handler {
    public:
        ReaderHandler(std::unique_ptr<Reader> &&reader, boost::dll::shared_library &library, std::shared_ptr<MessageChannel> channel)
            : reader(std::move(reader))
            , channel(std::move(channel))
            , library(library) {}

        void handle(std::iostream &stream) override {
            channel->push(reader->read(stream));
        }

        std::unique_ptr<Reader> reader;
        std::shared_ptr<MessageChannel> channel;
        boost::dll::shared_library library;
    };

    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //

    class Stream {
    public:
        struct entry {
            std::thread thread;
            std::unique_ptr<Node> node;
            boost::dll::shared_library library;
        };

        void join();

    private:
        std::vector<entry> entries;
    };

    void Stream::join() {
        for(auto &entry : this->entries) {
            entry.thread.join();
        }
    }

    class ConnectionImpl : public Connection, public std::enable_shared_from_this<ConnectionImpl> {
    public:

        using tcp = boost::asio::ip::tcp;

        ConnectionImpl(Context::Paths &paths, tcp::socket &socket)
                : stream(std::move(socket))
                , paths(paths)
                , channels {
                    std::make_shared<MessageChannel>(),
                    std::make_shared<MessageChannel>(),
                    std::make_shared<MessageChannel>()
                } {}

        void start();
        void process_input();
        void process_output();

        void build_stream(std::future<Config>, std::future<Header>);

        void initialize_readers(const Config &config);
        void initialize_writers(const Config &config);
        void initialize_stream(const Config &config, const Context &context);

        tcp::iostream stream;
        Gadgetron::Core::Context::Paths paths;

        struct {
            std::shared_ptr<MessageChannel> input, output, error;
        } channels;

        struct {
            std::promise<std::unique_ptr<Stream>> stream;

            std::promise<Config> config;
            std::promise<Header> header;

            std::promise<std::map<uint16_t, std::unique_ptr<Handler>>> readers;
            std::promise<std::vector<int>> writers;
        } promises;
    };


    void ConnectionImpl::start() {
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


    void ConnectionImpl::process_input() {

        GDEBUG_STREAM("Input thread running.");

        bool closed = false;

        auto reader_future = this->promises.readers.get_future();
        auto build_thread = std::thread(
                [&](auto config, auto header) { this->build_stream(std::move(config), std::move(header)); },
                this->promises.config.get_future(),
                this->promises.header.get_future()
        );
        build_thread.detach();

        std::map<uint16_t, std::unique_ptr<Handler>> handlers;
        handlers[FILENAME] = std::make_unique<ConfigFileHandler>(this->promises.config, paths);
        handlers[CONFIG]   = std::make_unique<ConfigHandler>(this->promises.config);
        handlers[HEADER]   = std::make_unique<HeaderHandler>(this->promises.header);
        handlers[CLOSE]    = std::make_unique<CloseHandler>(closed);
        handlers[QUERY]    = std::make_unique<QueryHandler>();

        while (!closed) {
            auto id = read_t<uint16_t>(stream);

            GDEBUG_STREAM("Handling message with id: " << id << std::endl);

            if (!handlers.count(id)) {
                handlers.merge(reader_future.get());
            }

            handlers.at(id)->handle(stream);
        }
    }

    void ConnectionImpl::process_output() {
        GDEBUG_STREAM("Output thread running.");
    }

    void ConnectionImpl::build_stream(std::future<Config> config_future, std::future<Header> header_future) {

        Config config = config_future.get();
        this->initialize_readers(config);
        this->initialize_writers(config);

        Context::Header header = header_future.get();
        Context context{header, paths};

        this->initialize_stream(config, context);
    }

    void ConnectionImpl::initialize_readers(const Config &config) {

        Builders::ReaderBuilder builder(config, paths);

        std::map<uint16_t, std::unique_ptr<Handler>> handlers;

        builder.process(
            [&](uint16_t port, std::unique_ptr<Reader> reader, boost::dll::shared_library lib) {
                handlers[port] = std::make_unique<ReaderHandler>(std::move(reader), lib, this->channels.input);
            }
        );

        this->promises.readers.set_value(std::move(handlers));
    }

    void ConnectionImpl::initialize_writers(const Config &config) {

        Builders::WriterBuilder builder(config, paths);

        // TODO: Writers
    }

    void ConnectionImpl::initialize_stream(const Config &config, const Context &context) {

        auto stream = std::make_unique<Stream>();

        Builders::StreamBuilder builder(config, context, channels.input, channels.output);

        builder.process(
            [&](std::unique_ptr<Node> node, boost::dll::shared_library lib) {
                GDEBUG_STREAM("Loaded node from shared library: " << lib.location() << std::endl);
            }
        );
    }
};

std::shared_ptr<Connection> Connection::create(Gadgetron::Core::Context::Paths &paths, tcp::socket &socket) {

    auto connection = std::make_shared<ConnectionImpl>(paths, socket);
    connection->start();

    return connection;
}
