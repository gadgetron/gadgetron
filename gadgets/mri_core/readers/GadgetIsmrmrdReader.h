#pragma  once
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/waveform.h>
#include <complex>
#include "Reader.h"
#include "Writer.h"
#include "mri_core_data.h"
#include "io/primitives.h"

namespace Gadgetron {


    /**
    Default implementation of GadgetMessageReader for IsmrmrdAcquisition messages
    */
class EXPORTGADGETSMRICORE GadgetIsmrmrdAcquisitionMessageReader : public Core::Reader {

    public:

        Core::Message read(std::istream& stream) final;


        uint16_t slot() final;
        ~GadgetIsmrmrdAcquisitionMessageReader() final = default;
    };

    // ------------------------------------------------------------------------------------------------------- //
    // ISMRMRD wave form reader/writer



class EXPORTGADGETSMRICORE GadgetIsmrmrdWaveformMessageReader : public Core::Reader {

    public:

        Core::Message read(std::istream& stream) final;
        uint16_t slot() final;
        ~GadgetIsmrmrdWaveformMessageReader() final = default;
    };
}
