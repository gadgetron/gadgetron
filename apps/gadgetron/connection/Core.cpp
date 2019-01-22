
#include "Core.h"

#include "ConfigConnection.h"
#include "Writers.h"

namespace {

    using namespace Gadgetron::Core;
    using namespace Gadgetron::Server::Connection;

    class RootErrorHandler : public ErrorHandler {
    public:

#if defined(NDEBUG)
        // When debugging, it is useful to have all exceptions bubble up to the
        // debugger. To enable this, we sabotage the error handler on debug builds.

        void handle(const std::string &location, std::function<void()> fn) override {
            try {
                fn();
            }
            catch (const ChannelClosed &e) {
                // Ignored.
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
            try {
                fn();
            }
            catch (const ChannelClosed &e) {
                // Ignored.
            }
        }
#endif

        void send_errors(std::iostream &stream) {

            Writers::TextWriter writer{};

            errors.close();
            InputChannel &errors(this->errors);

            for (auto error : errors) {
                writer.write(stream, error);
            }
        }

    private:
        void push_error(const std::string &location, const std::string &message) {
            std::string error("[" + location + "] ERROR: " + message);
            GERROR_STREAM(error);
            errors.push(error);
        }

        MessageChannel errors{};
    };


    void send_close(std::iostream &stream) {
        uint16_t close = 4;
        stream.write(reinterpret_cast<char *>(&close), sizeof(close));
    }

}


namespace Gadgetron::Server::Connection {

    void handle_connection(std::unique_ptr<std::iostream> stream, Gadgetron::Core::Context::Paths paths) {

        RootErrorHandler error_handler{};

        error_handler.handle("Connection Main Thread", [&]() {
            ConfigConnection::process(*stream, paths, error_handler);
        });

        error_handler.send_errors(*stream);
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