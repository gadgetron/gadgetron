#pragma once

#include "Writer.h"

#include "gadgetron_mricore_export.h"
#include <ismrmrd/waveform.h>
#include "hoNDArray.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GadgetIsmrmrdWaveformMessageWriter
            : public Core::TypedWriter<ISMRMRD::WaveformHeader, hoNDArray<uint32_t>> {

    public:
    protected:
        void serialize(std::ostream &stream, std::unique_ptr<ISMRMRD::WaveformHeader> header,
                       std::unique_ptr<hoNDArray<uint32_t>> array) override;
    };


    class EXPORTGADGETSMRICORE GadgetIsmrmrdAcquisitionMessageWriter
: public Core::TypedWriter<ISMRMRD::AcquisitionHeader, boost::optional<hoNDArray<float>>, hoNDArray<std::complex<float>>> {

    protected:
        void serialize(std::ostream &stream, std::unique_ptr<ISMRMRD::AcquisitionHeader> header ,std::unique_ptr<boost::optional<hoNDArray<float>>> trajectory , std::unique_ptr<hoNDArray<std::complex<float>>> data) override;
    };

}
