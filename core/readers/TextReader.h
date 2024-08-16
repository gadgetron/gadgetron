#pragma once

#include "Reader.h"

namespace Gadgetron::Core::Readers
{
    class TextReader : public Gadgetron::Core::Reader
    {
    public:
        Message read(std::istream &stream) override;
        uint16_t slot() override;
    };
}

