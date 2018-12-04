#ifndef GADGETRON_ACQUISITIONREADER_H
#define GADGETRON_ACQUISITIONREADER_H

#include "Reader.h"

namespace Gadgetron::Core::Readers {

    class AcquisitionReader : public Gadgetron::Core::Reader {
    public:
        virtual std::unique_ptr<Message> read(std::istream &stream) override;
        virtual uint16_t port() override;
    };
}

#endif //GADGETRON_ACQUISITIONREADER_H
