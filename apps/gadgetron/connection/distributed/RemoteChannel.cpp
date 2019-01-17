//
// Created by dchansen on 1/15/19.
//

#include <io/primitives.h>
#include "RemoteChannel.h"

#include <boost/asio.hpp>
#include "../Config.h"
#include "io/primitives.h"

using boost::asio::ip::tcp;
using namespace Gadgetron::Core;

#include "MessageID.h"

namespace {

    class MessageHandler {

    };

    void send_config_file(std::ostream &stream, const std::string &xml_config) {
        IO::write(stream, uint16_t(2));
        IO::write_string_to_stream<uint32_t>(stream, xml_config);
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

    if (!closed) IO::write(*stream, CLOSE);
    closed = true;
}

void Gadgetron::Server::Distributed::RemoteChannel::push_message(std::unique_ptr<Gadgetron::Core::Message> &&message) {

    std::lock_guard guard(closed_mutex);
    if (closed) throw Core::ChannelClosed();

    auto writer = std::find_if(writers.begin(), writers.end(),
                               [&](auto &writer) { return writer->accepts(*message); }
    );

    if (writer == writers.end())
        throw std::runtime_error("Tried to push message with no corresponding Writer to RemoteChannel");

    (*writer)->write(*stream, *message);

}

Gadgetron::Server::Distributed::RemoteChannel::RemoteChannel(const Address &address, const std::string &xml_config,
                                                             const std::map<uint16_t, std::shared_ptr<Gadgetron::Core::Reader>> &readers,
                                                             const std::vector<std::shared_ptr<Gadgetron::Core::Writer>> &writers)
        : readers(readers), writers(writers) {


    stream = std::make_unique<tcp::iostream>(address.ip, address.port);
    send_config_file(*stream, xml_config);


}

std::unique_ptr<Gadgetron::Core::Message> Gadgetron::Server::Distributed::RemoteChannel::pop() {
    auto id = Core::IO::read<uint16_t>(*stream);
    if (!is_writable_message(id)) throw error_readers.at(id)(*stream);

    return readers.at(id)->read(*stream);

}
