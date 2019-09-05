#pragma once

#include "Reader.h"
namespace Gadgetron::Core::Readers {

    class WaveformReader : public Gadgetron::Core::Reader {
    public:
        virtual Message read(std::istream &stream) override;
        virtual uint16_t slot() override;
    };
};
