#pragma  once
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "url_encode.h"
#include "gadgetron_mricore_export.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/waveform.h>
#include <complex>
#include "Reader.h"
#include "Writer.h"
#include "mri_core_data.h"
#include "io/readers.h"
#include "io/writers.h"





namespace Gadgetron {


    /**
    Default implementation of GadgetMessageReader for IsmrmrdAcquisition messages
    */
class EXPORTGADGETSMRICORE GadgetIsmrmrdAcquisitionMessageReader : public Core::Reader {

    public:

        virtual std::unique_ptr<Core::Message> read(std::istream& stream) override final;


        virtual uint16_t slot() override final;
        virtual ~GadgetIsmrmrdAcquisitionMessageReader(){};
    };

    // ------------------------------------------------------------------------------------------------------- //
    // ISMRMRD wave form reader/writer



class EXPORTGADGETSMRICORE GadgetIsmrmrdWaveformMessageReader : public Core::Reader {

    public:

        virtual std::unique_ptr<Core::Message> read(std::istream& stream) override final;
        virtual uint16_t slot() override final;
        virtual ~GadgetIsmrmrdWaveformMessageReader(){};
    };
}
