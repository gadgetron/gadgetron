
#include "RemoteChannel.h"
#include <io/primitives.h>

#include "../Config.h"
#include "io/primitives.h"
#include <boost/asio.hpp>

using boost::asio::ip::tcp;
using namespace Gadgetron::Core;

#include "MessageID.h"

namespace {

    void send_config_file(std::ostream& stream, const std::string& xml_config) {
        IO::write(stream, CONFIG);
        IO::write_string_to_stream<uint32_t>(stream, xml_config);
    }

    void send_header(std::ostream& stream, const ISMRMRD::IsmrmrdHeader& header) {

        std::stringstream sstream;
        ISMRMRD::serialize(header, sstream);
        IO::write(stream, HEADER);
        IO::write_string_to_stream<uint32_t>(stream, sstream.str());
    }

    bool is_writable_message(uint16_t message_id) {
        switch (message_id) {
        case FILENAME:
        case CONFIG:
        case HEADER:
        case CLOSE:
        case TEXT:
        case QUERY:
        case RESPONSE:
        case ERROR:
            return false;
        default:
            return true;
        }
    }

}

void Gadgetron::Server::Distributed::RemoteChannel::close() {

    std::lock_guard guard(closed_mutex);

    if (!closed_input)
        IO::write(*stream, CLOSE);
    closed_input = true;
}

void Gadgetron::Server::Distributed::RemoteChannel::push_message(Gadgetron::Core::Message message) {

    std::lock_guard guard(closed_mutex);
    if (closed_input)
        throw Core::ChannelClosed();

    try {
        auto writer
            = std::find_if(writers.begin(), writers.end(), [&](auto& writer) { return writer->accepts(message); });

        if (writer == writers.end())
            throw std::runtime_error("Tried to push message with no corresponding Writer to RemoteChannel");

        (*writer)->write(*stream, std::move(message));
    } catch (const std::ios_base::failure& failure) {
        throw RemoteError{ address, { failure.what() } };
    }
}

Gadgetron::Server::Distributed::RemoteChannel::RemoteChannel(const Address& address, const std::string& xml_config,
    const ISMRMRD::IsmrmrdHeader& header, const std::map<uint16_t, std::unique_ptr<Gadgetron::Core::Reader>>& readers,
    const std::vector<std::unique_ptr<Gadgetron::Core::Writer>>& writers)
    : readers(readers), writers(writers), address(address) {

    info_handlers = { { CLOSE, [this](auto& connection_stream) { this->handle_close(); } },
        { ERROR, [this](auto& connection_stream) {
             this->save_error(IO::read_string_from_stream<uint64_t>(connection_stream));
         } } };

    stream = std::make_unique<tcp::iostream>(tcp::v6(), address.ip, address.port);

    stream->exceptions(std::istream::failbit | std::istream::badbit | std::istream::eofbit);
    send_config_file(*stream, xml_config);
    send_header(*stream, header);
}

Gadgetron::Core::Message Gadgetron::Server::Distributed::RemoteChannel::pop() {
    {
        std::lock_guard guard(closed_mutex);
        if (closed_output)
            throw ChannelClosed();
    }
    auto id = Core::IO::read<uint16_t>(*stream);
    while (!is_writable_message(id)) {
        info_handlers.at(id)(*stream);
        id = Core::IO::read<uint16_t>(*stream);
    }
    return readers.at(id)->read(*stream);
}

void Gadgetron::Server::Distributed::RemoteChannel::handle_close() {
    std::lock_guard guard(closed_mutex);
    closed_input  = true;
    closed_output = true;
    if (error_messages.empty())
        throw ChannelClosed();
    throw RemoteError(address, error_messages);
}

void Gadgetron::Server::Distributed::RemoteChannel::save_error(const std::string& error_message) {
    error_messages.push_back(error_message);
}

namespace {
    std::string make_error_message(
        const Gadgetron::Server::Distributed::Address& address, const std::vector<std::string>& errors) {
        std::stringstream error_maker;
        error_maker << "Error received from " << address.ip << ":" << address.port << std::endl;
        error_maker << "Errors received: " << std::endl;
        for (auto& error : errors) {
            error_maker << error << std::endl;
        }
        return error_maker.str();
    }
}

Gadgetron::Server::Distributed::RemoteError::RemoteError(const Address& address, const std::vector<std::string>& errors)
    : std::runtime_error(make_error_message(address, errors)) {}
