
#include "io/primitives.h"
#include "mri_core_data.h"

#include "AcquisitionReader.h"

namespace Gadgetron::Core::Readers {


    Core::Message AcquisitionReader::read(std::istream &stream) {

        using namespace Core;
        using namespace std::literals;


        auto header = ISMRMRD::AcquisitionHeader{};

        IO::read(stream, header);

        optional<hoNDArray<float>> trajectory = boost::none;
        if (header.trajectory_dimensions) {
            trajectory = hoNDArray<float>(header.trajectory_dimensions,
                                          header.number_of_samples);
            IO::read(stream, *trajectory);
        }

        auto data = hoNDArray<std::complex<float>>(header.number_of_samples,
                                                   header.active_channels);
        IO::read(stream, data);

        return Core::Message(header, trajectory, data);
    }

    uint16_t AcquisitionReader::slot() {
        return 1008;
    }

    GADGETRON_READER_EXPORT(AcquisitionReader)
}


