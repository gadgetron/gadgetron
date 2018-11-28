#include <io/readers.h>
#include <boost/dll.hpp>

#include "hoNDArray.h"
#include "mri_core_data.h"

#include "AcquisitionReader.h"


namespace Gadgetron::Core::Readers {

     std::unique_ptr<Core::Message> AcquisitionReader::read(std::istream& stream) {

        using namespace Core;
        using namespace std::literals;

        auto acquisition = std::make_unique<Acquisition>();
        auto& header = acquisition->header;

        IO::read(stream, header);

        if (header.trajectory_dimensions) {
            acquisition->trajectory = hoNDArray<float>( header.trajectory_dimensions, header.number_of_samples);
            IO::read(stream,*acquisition->trajectory);
        }

        acquisition->data = hoNDArray<std::complex<float>>(header.number_of_samples, header.active_channels);
        IO::read(stream, acquisition->data);

        return std::unique_ptr<Message>(new TypedMessage<Acquisition>(std::move(acquisition)));
    }

    uint16_t AcquisitionReader::port() {
        return 1008;
    }

    std::shared_ptr<AcquisitionReader> reader_factory() {
        return std::make_shared<AcquisitionReader>();
    }
}

BOOST_DLL_ALIAS(
    Gadgetron::Core::Readers::reader_factory,
    reader_factory
);
