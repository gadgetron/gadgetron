#ifndef GADGETRON_ACQUISITIONREADER_H
#define GADGETRON_ACQUISITIONREADER_H

#include "Reader.h"

namespace Gadgetron::Core::Readers {

    class AcquisitionReader : public Gadgetron::Core::Reader {
    public:
        Message read(std::istream &stream) override;
        uint16_t slot() override;
    };
}

#endif //GADGETRON_ACQUISITIONREADER_H
