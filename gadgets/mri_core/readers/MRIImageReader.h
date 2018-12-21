#ifndef GADGETISMRMRDREADWRITE_H
#define GADGETISMRMRDREADWRITE_H

#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>

#include <complex>
#include "Reader.h"
#include "Message.h"

namespace Gadgetron{

    /**
    Default implementation of GadgetMessageReader for Ismrmrd Images with attribues messages
    */
class EXPORTGADGETSMRICORE MRIImageReader : public Core::Reader
    {
    public:

        std::unique_ptr<Core::Message> read(std::istream& stream) override;
        uint16_t slot() override;

    };

}
#endif //GADGETISMRMRDREADWRITE_H
