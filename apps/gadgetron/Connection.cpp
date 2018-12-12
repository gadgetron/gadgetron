
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

#include "writers/ResponseWriter.h"
#include "Response.h"
#include "Builders.h"
#include "Reader.h"
#include "Writer.h"
#include "Server.h"
#include "Config.h"
#include "Query.h"

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
        QUERY       = 6,
        RESPONSE    = 7
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

    template<class T>
    std::string read_string_from_stream(std::istream &stream) {

        auto n = read_t<T>(stream);
        auto buffer = std::make_unique<char[]>(n);

        stream.read(buffer.get(), n);

        return std::string(buffer.get());
    }

    class Handler {
    public:
        virtual void handle(std::istream &stream) = 0;
        virtual ~Handler() = default;
    };

    class ConfigFileHandler : public Handler {
    public:
        ConfigFileHandler(
                std::promise<Config> &config_promise,
                std::promise<std::stringstream> &raw_config_promise,
                const Context::Paths &paths)
                : promises {config_promise, raw_config_promise}
                , paths(paths) {}

        void handle(std::istream &stream) override {
            boost::filesystem::path filename = paths.gadgetron_home / GADGETRON_CONFIG_PATH / read_filename_from_stream(stream);

            GDEBUG_STREAM("Reading config file: " << filename << std::endl);

            std::ifstream file_stream(filename.string());

            std::stringstream config_stream;
            config_stream << file_stream.rdbuf();

            promises.config.set_value(parse_config(config_stream));
            promises.string.set_value(std::move(config_stream));
        }


    private:
        struct {
            std::promise<Config> &config;
            std::promise<std::stringstream> &string;
        } promises;

        const Context::Paths &paths;
    };

    class ConfigHandler : public Handler {
    public:
        explicit ConfigHandler(
                std::promise<Config> &config_promise,
                std::promise<std::stringstream> &raw_config_promise
        ) : promises {config_promise, raw_config_promise} {}

        void handle(std::istream &stream) override {
            std::stringstream config_stream(read_string_from_stream<uint32_t>(stream));

            promises.config.set_value(parse_config(config_stream));
            promises.string.set_value(std::move(config_stream));
        }

    private:
        struct {
            std::promise<Config> &config;
            std::promise<std::stringstream> &string;
        } promises;
    };

    class HeaderHandler : public Handler {
    public:
        explicit HeaderHandler(std::promise<Header> &header_promise) : promise(header_promise) {}

        void handle(std::istream &stream) override {
            std::string raw_header(read_string_from_stream<uint32_t>(stream));

            ISMRMRD::IsmrmrdHeader header;
            ISMRMRD::deserialize(raw_header.c_str(), header);

            promise.set_value(header);
        }

    private:
        std::promise<Header> &promise;
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

    class QueryHandler : public Handler {
    public:
        explicit QueryHandler(
                std::future<std::stringstream> &&raw_config_future,
                std::shared_ptr<MessageChannel> output
        ) : channel(std::move(output)) {
            handlers.emplace_back(std::make_unique<Query::ISMRMRDHandler>());
            handlers.emplace_back(std::make_unique<Query::GadgetronHandler>());
            handlers.emplace_back(std::make_unique<Query::ConfigHandler>(std::move(raw_config_future)));
        }



        void handle(std::istream &stream) override {

            auto reserved = read_t<uint64_t>(stream);
            auto corr_id  = read_t<uint64_t>(stream);
            auto query    = read_string_from_stream<uint64_t>(stream);

            if (reserved) {
                throw std::runtime_error("Unsupported value in reserved bytes.");
            }

            auto response = get_response(query);

            channel->push(std::make_unique<Response>(corr_id, response));
        }

        std::string get_response(const std::string &query) {

            GDEBUG_STREAM("Processing query: " << query);

            for (auto &handler : handlers) {
                if (handler->accepts(query)) {
                    return handler->handle(query);
                }
            }

            // Response should be "Unable to respond."

            throw std::runtime_error("No query handler found for query: " + query);
        }

    private:
        std::vector<std::unique_ptr<Query::Handler>> handlers;
        std::shared_ptr<OutputChannel> channel;
    };

    class ReaderHandler : public Handler {
    public:
        ReaderHandler(std::unique_ptr<Reader> &&reader, std::shared_ptr<MessageChannel> channel)
                : reader(std::move(reader)), channel(std::move(channel)) {}

        void handle(std::istream &stream) override {
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

        ConnectionImpl(Context::Paths &paths_in, std::unique_ptr<tcp::iostream>& stream_in)
                : stream(std::move(stream_in))
                , paths(paths_in)
                , builder(paths_in)
                , channels {
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

        std::promise<std::stringstream> raw_config_promise;

        auto reader_future = this->promises.readers.get_future();
        auto stream_thread = std::thread(
                [&](auto config, auto header) { this->start_stream(std::move(config), std::move(header)); },
                this->promises.config.get_future(),
                this->promises.header.get_future()
        );

        std::map<uint16_t, std::unique_ptr<Handler>> handlers;
        handlers[FILENAME] = std::make_unique<ConfigFileHandler>(this->promises.config, raw_config_promise, paths);
        handlers[CONFIG]   = std::make_unique<ConfigHandler>(this->promises.config, raw_config_promise);
        handlers[HEADER]   = std::make_unique<HeaderHandler>(this->promises.header);
        handlers[CLOSE]    = std::make_unique<CloseHandler>(closed);
        handlers[QUERY]    = std::make_unique<QueryHandler>(raw_config_promise.get_future(), channels.output);

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

//        auto writers = std::vector<std::unique_ptr<Writer>>();
//        writers.push_back(std::make_unique<Writers::ResponseWriter>());

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
