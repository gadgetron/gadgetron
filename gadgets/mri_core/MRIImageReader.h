#ifndef GADGETISMRMRDREADWRITE_H
#define GADGETISMRMRDREADWRITE_H

#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "GadgetMessageInterface.h"
#include "hoNDArray.h"
#include "url_encode.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>

#include <ace/SOCK_Stream.h>
#include <ace/Task.h>
#include <complex>

namespace Gadgetron{

    /**
    Default implementation of GadgetMessageReader for Ismrmrd Images with attribues messages
    */
    class EXPORTGADGETSMRICORE MRIImageReader : public GadgetMessageReader
    {
    public:
        GADGETRON_READER_DECLARE(MRIImageReader);
        virtual ACE_Message_Block* read(ACE_SOCK_Stream* stream);
    };

}
#endif //GADGETISMRMRDREADWRITE_H
