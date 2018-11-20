#ifndef GADGETRON_STREAM_H
#define GADGETRON_STREAM_H

#include <memory>
#include "Message.h"
namespace Gadgetron::Core {

    class Stream {

    };

    class Reader {
    public:
        virtual std::unique_ptr<Message> read(std::istream &stream);
    };

    class Writer {
    public:
        virtual void write(std::ostream &stream, std::unique_ptr<Message>&& message);
    };
}



#endif //GADGETRON_STREAM_H
