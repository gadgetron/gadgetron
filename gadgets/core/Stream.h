#pragma once
#include <memory>
#include <typeindex>
#include <boost/dll.hpp>

#include "Message.h"
#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Core {

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

#define GADGETRON_READER_EXPORT(ReaderClass)                        \
std::unique_ptr<Reader> reader_factory_##ReaderClass() {            \
    return std::make_unique<ReaderClass>();                         \
}                                                                   \
                                                                    \
BOOST_DLL_ALIAS(                                                    \
        reader_factory_##ReaderClass,                               \
        reader_factory_export_##ReaderClass                         \
)                                                                   \

#define GADGETRON_WRITER_EXPORT(WriterClass)

#define GADGETRON_NODE_EXPORT(NodeClass)                            \


#define GADGETRON_GADGET_EXPORT(GadgetClass)



