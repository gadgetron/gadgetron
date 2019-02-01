
#include "Core.h"

#include "ConfigConnection.h"
#include "Writers.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Server::Connection;

    class ErrorSender : public ErrorReporter {
    public:


        void operator()(const std::string &location, const std::string &message) override {
            std::string error("[" + location + "] ERROR: " + message);
            GERROR_STREAM(error);
            {
                std::lock_guard guard(error_lock);
                errors.push_back(error);
            }
        }


        void send_errors(std::iostream &stream) {

            Writers::TextWriter writer{};

            for (auto &error : errors) {
                writer.serialize(stream, error);
            }
        }

    private:

        std::mutex error_lock;
        std::vector<std::string> errors;

    };


    void send_close(std::iostream &stream) {
        uint16_t close = 4;
        stream.write(reinterpret_cast<char *>(&close), sizeof(close));
    }

}


namespace Gadgetron::Server::Connection {

    void handle_connection(std::unique_ptr<std::iostream> stream, Gadgetron::Core::Context::Paths paths) {

        stream->exceptions(std::istream::failbit | std::istream::badbit | std::istream::eofbit);
        ErrorSender sender;

        ErrorHandler error_handler(sender,"Connection Main Thread");

        error_handler.handle([&]() {
            ConfigConnection::process(*stream, paths, error_handler);
        });

        sender.send_errors(*stream);
        send_close(*stream);
    }

    std::vector<std::unique_ptr<Core::Writer>> default_writers() {
        std::vector<std::unique_ptr<Writer>> writers{};

        writers.emplace_back(std::make_unique<Writers::TextWriter>());
        // writers.emplace_back(std::make_unique<Writers::ErrorWriter>());
        writers.emplace_back(std::make_unique<Writers::ResponseWriter>());

        return std::move(writers);
    }
}