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
        virtual Message read(std::istream &stream) = 0;
        virtual uint16_t slot() = 0;

        virtual ~Reader() = default;
    };
}

#define GADGETRON_READER_EXPORT(ReaderClass)                        \
std::unique_ptr<Gadgetron::Core::Reader> reader_factory_##ReaderClass() {            \
    return std::make_unique<ReaderClass>();                         \
}                                                                   \
                                                                    \
BOOST_DLL_ALIAS(                                                    \
        reader_factory_##ReaderClass,                               \
        reader_factory_export_##ReaderClass                         \
)                                                                   \

