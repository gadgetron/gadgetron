#ifndef GADGETRON_WAVEFORMREADER_H
#define GADGETRON_WAVEFORMREADER_H

#include "Reader.h"

namespace Gadgetron::Core::Readers {

    class WaveformReader : public Gadgetron::Core::Reader {
    public:
        virtual std::unique_ptr<Message> read(std::istream &stream) override;
        virtual uint16_t slot() override;
    };
};


#endif //GADGETRON_WAVEFORMREADER_H
