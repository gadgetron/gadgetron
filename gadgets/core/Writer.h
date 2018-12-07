#pragma once
#include <memory>
#include <typeindex>
#include <boost/dll.hpp>

#include "Message.h"
#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Core {

    class Writer {
    public:
        virtual ~Writer() = default;
        virtual bool accepts(const Message &) = 0;
        virtual void write(std::ostream &stream, std::unique_ptr<Message> &&message) = 0;
    };

    template<class ...ARGS>
    class TypedWriter : public Writer {
    public:
        ~TypedWriter() override = default;
        bool accepts(const Message &) override;
        void write(std::ostream &stream, std::unique_ptr<Message> &&message) override;

    protected:
        virtual void serialize(std::ostream &stream, std::unique_ptr<ARGS> &&...) = 0;
    };

}

namespace Gadgetron::Core {

    template<class ...ARGS>
    bool TypedWriter<ARGS...>::accepts(const Message &message) {
        // TODO: Then magic happens.

        return false;
    }

    template<class ...ARGS>
    void TypedWriter<ARGS...>::write(std::ostream &stream, std::unique_ptr<Message> &&message) {

    }
}


#include "Writer.hpp"

#define GADGETRON_WRITER_EXPORT(WriterClass)                /
