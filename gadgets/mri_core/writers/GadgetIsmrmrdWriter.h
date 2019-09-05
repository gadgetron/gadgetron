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
        void serialize(std::ostream &stream, const ISMRMRD::WaveformHeader &header,
                       const hoNDArray<uint32_t> &array) override;
    };


    class EXPORTGADGETSMRICORE GadgetIsmrmrdAcquisitionMessageWriter
            : public Core::TypedWriter<ISMRMRD::AcquisitionHeader, boost::optional<hoNDArray<float>>, hoNDArray<std::complex<float>>> {

    protected:
        void serialize(std::ostream &stream, const ISMRMRD::AcquisitionHeader &header,
                       const boost::optional<hoNDArray<float>> &trajectory,
                       const hoNDArray<std::complex<float>> &data) override;
    };

}
