#include "ImageAccumulator.h"

#include <chrono>
#include <ismrmrd/ismrmrd.h>
#include <boost/range/algorithm/copy.hpp>

#include "Unmixing.h"
#include "common/AcquisitionBuffer.h"

#include "Node.h"
#include "hoNDArray.h"
#include "hoNDFFT.h"

#include "log.h"
#include <range/v3/action.hpp>

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;


    Grappa::Image create_reconstruction_job(
            const AnnotatedAcquisition &first,
            const AnnotatedAcquisition &last,
            AcquisitionBuffer &buffer
    ) {
        auto &first_header = std::get<ISMRMRD::AcquisitionHeader>(first);
        auto &last_header = std::get<ISMRMRD::AcquisitionHeader>(last);
        auto slice = last_header.idx.slice;

        Grappa::Image image{};

        image.data = buffer.take(slice);
        image.meta.slice = slice;
        image.meta.time_stamp = first_header.acquisition_time_stamp;

        boost::copy(last_header.position, image.meta.position.begin());
        boost::copy(last_header.read_dir, image.meta.read_dir.begin());
        boost::copy(last_header.phase_dir, image.meta.phase_dir.begin());
        boost::copy(last_header.slice_dir, image.meta.slice_dir.begin());
        boost::copy(last_header.patient_table_position, image.meta.table_pos.begin());

        hoNDFFT<float>::instance()->ifft2c(image.data);

        return std::move(image);
    }
}

namespace Gadgetron::Grappa {

    ImageAccumulator::ImageAccumulator(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : ChannelGadget<Slice>(context,props), context(context) {}

    void ImageAccumulator::process(InputChannel<Slice> &in, OutputChannel &out) {

        AcquisitionBuffer buffer{context};

        for (auto slice : in) {

            slice = std::move(slice) | ranges::actions::remove_if([](auto& acq){return std::get<ISMRMRD::AcquisitionHeader>(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION);});

            buffer.add(slice);
            out.push(create_reconstruction_job(slice.front(), slice.back(), buffer));
        }
    }

    GADGETRON_GADGET_EXPORT(ImageAccumulator);
}
