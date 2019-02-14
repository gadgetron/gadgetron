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

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;


    void emit_reconstruction_job(
            const AnnotatedAcquisition &acquisition,
            AcquisitionBuffer &buffer,
            OutputChannel &output
    ) {
        auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        auto slice = header.idx.slice;

        Grappa::Image image{};

        image.data = buffer.take(slice);
        image.meta.slice = slice;

        boost::copy(header.position, image.meta.position.begin());
        boost::copy(header.read_dir, image.meta.read_dir.begin());
        boost::copy(header.phase_dir, image.meta.phase_dir.begin());
        boost::copy(header.slice_dir, image.meta.slice_dir.begin());
        boost::copy(header.patient_table_position, image.meta.table_pos.begin());

        hoNDFFT<float>::instance()->ifft3c(image.data);

        output.push(std::move(image));
    }
}

namespace Gadgetron::Grappa {

    ImageAccumulator::ImageAccumulator(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode<Slice>(props), context(context) {}

    void ImageAccumulator::process(TypedInputChannel<Slice> &in, OutputChannel &out) {

        AcquisitionBuffer buffer{context};

        for (const auto &slice : in) {
            buffer.add(slice);
            emit_reconstruction_job(slice.back(), buffer, out);
        }
    }

    GADGETRON_GADGET_EXPORT(ImageAccumulator);
}
