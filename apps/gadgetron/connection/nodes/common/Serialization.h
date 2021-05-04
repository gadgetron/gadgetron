#pragma once

#include <map>
#include <memory>
#include <functional>

#include "Reader.h"
#include "Writer.h"

namespace Gadgetron::Server::Connection::Nodes {

    class Serialization {
    public:
        using Readers = std::map<uint16_t, std::unique_ptr<Core::Reader>>;
        using Writers = std::vector<std::unique_ptr<Core::Writer>>;
        Serialization(Readers readers, Writers writers);

        void close(std::iostream &stream) const;
        void write(std::iostream &stream, Core::Message message) const;
        Core::Message read(
                std::iostream &stream,
                std::function<void()> on_close,
                std::function<void(std::string message)> on_error
        ) const;
    private:
        const Readers readers;
        const Writers writers;
    };
}