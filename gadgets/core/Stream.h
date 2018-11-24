#ifndef GADGETRON_STREAM_H
#define GADGETRON_STREAM_H

#include <memory>
#include "Message.h"
#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Core {

    class Stream : public Node {
    public:
        Stream(const std::string &stream_name, std::shared_ptr<InputChannel<Message>> channel) : name(stream_name),
                                                                                                 output(channel) {}

        virtual std::shared_ptr<InputChannel<Message>> output_channel() override { return output;}
    protected:
        std::string name;
        std::shared_ptr<InputChannel<Message>> output;


    };

    class Reader {
    public:
        virtual std::unique_ptr<Message> read(std::istream &stream) = 0;

        virtual uint16_t port() = 0;
    };

    class Writer {
    public:
        virtual void write(std::ostream &stream, std::unique_ptr<Message> &&message) = 0;
    };
}


#endif //GADGETRON_STREAM_H
