#pragma once

#include "Reader.h"
#include "MessageID.h"

namespace Gadgetron::Core::Readers {
    class BufferReader : public Gadgetron::Core::Reader {
    public:
        Message read(std::istream &stream) override;
        uint16_t slot() override;
    };
}


