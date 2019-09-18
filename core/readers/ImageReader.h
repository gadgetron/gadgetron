//
// Created by dchansen on 9/13/19.
//

#pragma once
#include "Reader.h"

namespace Gadgetron::Core::Readers {
    class ImageReader : public Core::Reader {
    public:
        Message read(std::istream& stream) override;
        uint16_t slot() override;
    };
}
