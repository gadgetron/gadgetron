#pragma once
#include <memory>
#include <typeindex>
#include "Message.h"
#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Core {

    class Stream {
    public:
        Stream(const std::string &stream_name, const std::shared_ptr<InputChannel<Message>>& channel) : name(stream_name),
                                                                                                 output(channel) {}

        virtual std::shared_ptr<InputChannel<Message>> output_channel() { return output;}
    protected:
        std::string name;
        std::shared_ptr<InputChannel<Message>> output;


    };

    class Reader {
    public:
        virtual std::unique_ptr<Message> read(std::istream &stream) = 0;
        virtual uint16_t port() = 0;
        virtual ~Reader() {};
    };

    class Writer {
    public:
        virtual void write(std::ostream &stream, std::unique_ptr<Message> &&message) = 0;
        virtual std::vector<std::type_index> supported_types() const = 0;
        virtual ~Writer() {};
    };
}


