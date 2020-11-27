#pragma once
#include "Reader.h"

namespace Gadgetron::Core::Readers {
    class IsmrmrdImageArrayReader : public Reader {
    public:
        Message read(std::istream& stream) override;
        uint16_t slot() override;
    };
}
