
#include <memory>

#include "connection/VoidConnection.h"
#include "connection/ConfigConnection.h"
#include "connection/HeaderConnection.h"
#include "connection/StreamConnection.h"
#include "connection/Writers.h"

#include "Connection.h"


#include "readers/Primitives.h"
#include "log.h"

using namespace boost::asio;

using namespace Gadgetron::Core;
using namespace Gadgetron::Core::Readers;

using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Writers;

namespace {

    template<class T>
    struct RAIICloser {
        std::shared_ptr<T> t;
        ~RAIICloser() { t->close(); }
    };

    class ErrorChannel : public ErrorHandler {
    public:

#ifdef NDEBUG
        // When debugging, it is useful to have all exceptions bubble up to the
        // debugger. To enable this, we sabotage the error handler on debug builds.

        void handle(const std::string &location, std::function<void()> fn) override {
            try {
                fn();
            }
            catch (const std::exception &e) {
                push_error(location, e.what());
            }
            catch (const std::string &s) {
                push_error(location, s);
            }
            catch (...) {
                push_error(location, "Unknown error.");
            }
        }
#else
        void handle(const std::string &, std::function<void()> fn) override {
            fn();
        }
#endif

        void send_errors(std::iostream &stream) {

            TextWriter writer{};

            errors.close();
            InputChannel<Message> &errors(this->errors);

            for (auto error : errors) {
                writer.write(stream, std::move(error));
            }
        }

    private:
        void push_error(const std::string &location, const std::string &message) {
            errors.push(std::make_unique<std::string>("[" + location + "] ERROR: " + message));
        }

        MessageChannel errors{};
    };


    void send_close(std::iostream &stream) {
        uint16_t close = 4;
        stream.write(reinterpret_cast<char *>(&close), sizeof(close));
    }

    void handle_connection(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {

        ErrorChannel error_handler{};

        error_handler.handle("Connection Main Thread", [&]() {
            ConfigConnection::process(*stream, paths, error_handler);
        });

        error_handler.send_errors(*stream);
        send_close(*stream);
    }
}

namespace Gadgetron::Server::Connection {

    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {
        auto thread = std::thread(handle_connection, paths, std::move(stream));
        thread.detach();
    }

    Connection::Connection(std::iostream &stream) : stream(stream) {};

    void Connection::process_input() {

        RAIICloser<MessageChannel> closer{channels.input};

        bool closed = false;
        auto handlers = prepare_handlers([&]() { closed = true; });

        while (!closed) {
            auto id = read_t<uint16_t>(stream);

            GDEBUG_STREAM("Processing message with id: " << id);

            handlers.at(id)->handle(stream);
        }
    }

    void Connection::process_output() {

        RAIICloser<MessageChannel> closer{channels.output};

        auto writers = prepare_writers();

        writers.push_back(std::make_unique<TextWriter>());
        writers.push_back(std::make_unique<ResponseWriter>());

        InputChannel<Message> &output = *channels.output;
        for (auto message : output) {

            auto writer = std::find_if(writers.begin(), writers.end(),
                   [&](auto &writer) { return writer->accepts(*message); }
            );

            if (writer != writers.end()) {
                (*writer)->write(stream, std::move(message));
            }
        }
    }

    std::vector<std::unique_ptr<Writer>> Connection::prepare_writers() {
        return std::vector<std::unique_ptr<Writer>>();
    }
}

