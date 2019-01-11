

#include "Writers.h"

namespace Gadgetron::Server::Connection::Writers {

    void ResponseWriter::serialize(
            std::ostream &stream,
            const Gadgetron::Core::Response& response
    ) {
        uint16_t message_id = 7;
        uint64_t correlation_id = response.correlation_id;
        uint64_t response_length = response.response.size();

        stream.write(reinterpret_cast<char *>(&message_id), sizeof(message_id));
        stream.write(reinterpret_cast<char *>(&correlation_id), sizeof(correlation_id));
        stream.write(reinterpret_cast<char *>(&response_length), sizeof(response_length));
        stream.write(response.response.c_str(), response.response.size());
    }


    void TextWriter::serialize(
            std::ostream &stream,
            const std::string& message
    ) {
        uint16_t message_id = 5;

        auto length = uint32_t(message.size());

        stream.write(reinterpret_cast<char *>(&message_id), sizeof(message_id));
        stream.write(reinterpret_cast<char *>(&length), sizeof(length));
        stream.write(message.data(), length);
    }
}
