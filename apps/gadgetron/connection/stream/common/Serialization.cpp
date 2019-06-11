
#include <iomanip>

#include "Serialization.h"

#include "io/primitives.h"
#include "MessageID.h"

using namespace Gadgetron::Core;

namespace Gadgetron::Server::Connection::Stream {

    Serialization::Serialization(
            Readers readers,
            Writers writers
    ) : readers(std::move(readers)), writers(std::move(writers)) {}

    void Serialization::write(std::iostream &stream, Core::Message message) const {

        auto writer = std::find_if(
                writers.begin(),
                writers.end(),
                [&](auto& writer) { return writer->accepts(message); }
        );

        if (writer == writers.end())
            throw std::runtime_error("Could not find appropriate writer for message.");

        (*writer)->write(stream, std::move(message));
    }

    Core::Message Serialization::read(
            std::iostream &stream,
            std::function<void()> on_close,
            std::function<void(std::string message)> on_error
    ) const {

        auto id = IO::read<uint16_t>(stream);
        auto illegal_message = [&](auto &) {
            throw std::runtime_error("Received illegal message id from external peer: " + std::to_string(id));
        };

        std::map<uint16_t, std::function<void(std::iostream &)>> handlers{
                {FILENAME,  illegal_message},
                {CONFIG,    illegal_message},
                {HEADER,    illegal_message},
                {CLOSE,     [&](auto &) { on_close(); }},
                {TEXT,      illegal_message},
                {QUERY,     illegal_message},
                {RESPONSE,  illegal_message},
                {ERROR,     [&](auto &stream) { on_error(IO::read_string_from_stream<uint64_t>(stream)); }}
        };

        for (; handlers.count(id); id = IO::read<uint16_t>(stream)) handlers.at(id)(stream);

        if (!readers.count(id)) {
            GERROR_STREAM("Oh shit, someone gave us a retarded message id: " << id);

            auto bytes = std::make_unique<char[]>(10);
            stream.read(bytes.get(), 10);

            std::stringstream output;
            output << std::hex << std::setw(2) << std::setfill('0');

            for (int i = 0; i < 10; i++) {
                output << "0x" << int(uint8_t(bytes[i])) << " ";
            }

            GERROR_STREAM("Let's have a look at the next few bytes: " << output.str());
        }

        return readers.at(id)->read(stream);
    }

    void Serialization::close(std::iostream &stream) const {
        IO::write(stream, CLOSE);
    }
}