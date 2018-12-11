
#include <typeindex>
#include <iostream>
#include <sstream>
#include <string>
#include <boost/dll/shared_library.hpp>
#include <boost/dll.hpp>
#include <boost/range/algorithm/transform.hpp>

#include "log.h"
#include "gadgetron_config.h"

#include "Connection.h"

#include "Builders.h"
#include "Reader.h"
#include "Writer.h"
#include "Server.h"
#include "Config.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Server;

    using Header = Gadgetron::Core::Context::Header;

    enum message_id : uint16_t {
        FILENAME = 1,
        CONFIG = 2,
        HEADER = 3,
        CLOSE = 4,
        TEXT = 5,
        QUERY = 6
    };

    template<class T>
    void read_into(std::istream &stream, T &t) {
        stream.read(reinterpret_cast<char *>(&t), sizeof(t));
    }

    template<class T>
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

        virtual ~Handler() {};
    };

    class ConfigFileHandler : public Handler {
    public:
        ConfigFileHandler(std::promise<Gadgetron::Server::Config> &config_promise, const Context::Paths &paths)
                : promise(config_promise), paths(paths) {}

        void handle(std::iostream &stream) override {
            boost::filesystem::path filename =
                    paths.gadgetron_home / GADGETRON_CONFIG_PATH / read_filename_from_stream(stream);
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
        ConfigHandler(std::promise<Config> &config_promise) : promise(config_promise) {}

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
        ReaderHandler(std::unique_ptr<Reader> &&reader, std::shared_ptr<MessageChannel> channel)
                : reader(std::move(reader)), channel(std::move(channel)) {}

        void handle(std::iostream &stream) override {
            channel->push_message(reader->read(stream));
        }

        virtual ~ReaderHandler() {};
        std::unique_ptr<Reader> reader;
        std::shared_ptr<MessageChannel> channel;
    };

    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //



    class ConnectionImpl : public Connection, public std::enable_shared_from_this<ConnectionImpl> {
    public:

        using tcp = boost::asio::ip::tcp;

        ConnectionImpl(Context::Paths &paths_in, std::unique_ptr<tcp::iostream> &stream_in)
                : stream(std::move(stream_in)), paths(paths_in), builder(paths_in), channels{
                std::make_shared<MessageChannel>(),
                std::make_shared<MessageChannel>()
        } {}

        void start();

        void process_input();

        void process_output();

        void start_stream(std::future<Config>, std::future<Header>);

        void initialize_readers(const Config &config);

        void initialize_writers(const Config &config);

        std::unique_ptr<tcp::iostream> stream;
        const Gadgetron::Core::Context::Paths paths;
        Builder builder;

        struct {
            std::shared_ptr<MessageChannel> input, output;
        } channels;

        struct {
            std::promise<std::unique_ptr<Gadgetron::Core::Node>> stream;

            std::promise<Config> config;
            std::promise<Header> header;

            std::promise<std::map<uint16_t, std::unique_ptr<Handler>>> readers;
            std::promise<std::vector<std::unique_ptr<Writer>>> writers;
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
        auto stream_thread = std::thread(
                [&](auto config, auto header) { this->start_stream(std::move(config), std::move(header)); },
                this->promises.config.get_future(),
                this->promises.header.get_future()
        );

        std::map<uint16_t, std::unique_ptr<Handler>> handlers;
        handlers[FILENAME] = std::make_unique<ConfigFileHandler>(this->promises.config, paths);
        handlers[CONFIG] = std::make_unique<ConfigHandler>(this->promises.config);
        handlers[HEADER] = std::make_unique<HeaderHandler>(this->promises.header);
        handlers[CLOSE] = std::make_unique<CloseHandler>(closed);
        handlers[QUERY] = std::make_unique<QueryHandler>();

        while (!closed) {
            auto id = read_t<uint16_t>(*stream);

            GDEBUG_STREAM("Handling message with id: " << id << std::endl);

            if (!handlers.count(id)) {
                handlers.merge(reader_future.get());
            }

            handlers.at(id)->handle(*stream);
        }
        channels.input->close();
        stream_thread.join();
    }

    void ConnectionImpl::process_output() {
        GDEBUG_STREAM("Output thread running.");

        auto writer_future = this->promises.writers.get_future();
        auto writers = writer_future.get();

        std::shared_ptr<InputChannel<Message>> output = this->channels.output;

//        for (std::unique_ptr<Message> message : *output) {
        try {
            while (true) {
                auto message = output->pop();
                GDEBUG_STREAM("Ptr " << message.get() << std::endl);
                GDEBUG_STREAM("Writer got a: " << typeid(*message).name() << std::endl);
                auto writer = std::find_if(writers.begin(), writers.end(),
                                           [&](auto &writer) { return writer->accepts(*message); }
                );


                if (writer != writers.end()) {
                    (*writer)->write(*stream, std::move(message));
                }
            }
        } catch (ChannelClosedError err) {

        }
    }


    void ConnectionImpl::start_stream(std::future<Config> config_future, std::future<Header> header_future) {

        Config config = config_future.get();
        this->initialize_readers(config);
        this->initialize_writers(config);

        Context::Header header = header_future.get();
        Context context{header, paths};

        auto stream = builder.build_stream(config.stream, context);

        stream->process(channels.input, channels.output);
    }

    void ConnectionImpl::initialize_readers(const Config &config) {


        std::map<uint16_t, std::unique_ptr<Handler>> handlers;

        auto readers = builder.build_readers(config.readers);

        for (auto &reader_pair : readers) {
            handlers.emplace(reader_pair.first,
                             std::make_unique<ReaderHandler>(std::move(reader_pair.second), channels.input));
        }
        this->promises.readers.set_value(std::move(handlers));
    }

    void ConnectionImpl::initialize_writers(const Config &config) {
        this->promises.writers.set_value(builder.build_writers(config.writers));

    }


};

std::shared_ptr<Connection>
Connection::create(Gadgetron::Core::Context::Paths &paths, std::unique_ptr<tcp::iostream> &stream) {

    auto connection = std::make_shared<ConnectionImpl>(paths, stream);
    connection->start();

    return connection;
}
