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
        virtual void write(std::ostream &stream, std::unique_ptr<Message> &&message) = 0;
        virtual std::vector<std::type_index> supported_types() const = 0;
        virtual ~Writer() {};
    };
}

#define GADGETRON_WRITER_EXPORT(WriterClass)
