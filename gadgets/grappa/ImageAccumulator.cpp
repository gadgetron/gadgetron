#include "ImageAccumulator.h"

#include <chrono>
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
    using namespace Gadgetron::Grappa;


    Grappa::Image create_reconstruction_job(
            const AnnotatedAcquisition &first,
            const AnnotatedAcquisition &last,
            AcquisitionBuffer &buffer
    ) {
        auto &first_acq = std::get<mrd::Acquisition>(first);
        auto &last_acq = std::get<mrd::Acquisition>(last);
        auto slice = last_acq.head.idx.slice.value_or(0);

        Grappa::Image image{};

        image.data = buffer.take(slice);
        image.meta.slice = slice;
        image.meta.time_stamp = first_acq.head.acquisition_time_stamp.value_or(0);

        boost::copy(last_acq.head.position, image.meta.position.begin());
        boost::copy(last_acq.head.read_dir, image.meta.read_dir.begin());
        boost::copy(last_acq.head.phase_dir, image.meta.phase_dir.begin());
        boost::copy(last_acq.head.slice_dir, image.meta.slice_dir.begin());
        boost::copy(last_acq.head.patient_table_position, image.meta.table_pos.begin());

        hoNDFFT<float>::instance()->ifft2c(image.data);

        return std::move(image);
    }
}

namespace Gadgetron::Grappa {

    ImageAccumulator::ImageAccumulator(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : ChannelGadget<Slice>(context,props), context(context) {}

    void ImageAccumulator::process(Core::InputChannel<Slice> &in, Core::OutputChannel &out) {

        AcquisitionBuffer buffer{context};

        for (auto slice : in) {

            slice = std::move(slice) | ranges::actions::remove_if([](auto& acq){return std::get<mrd::Acquisition>(acq).head.flags.HasFlags(mrd::AcquisitionFlags::kIsParallelCalibration);});

            buffer.add(slice);
            out.push(create_reconstruction_job(slice.front(), slice.back(), buffer));
        }
    }

    GADGETRON_GADGET_EXPORT(ImageAccumulator);
}
