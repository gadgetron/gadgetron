#pragma once

#include <cpu/hoNDArray.h>

#include "Writer.h"

namespace Gadgetron::Core::Writers {

    class ImageWriter : public Writer {
    public:
        bool accepts(const Message &) override;
        void write(std::ostream &stream, Message message) override;
    };
}

