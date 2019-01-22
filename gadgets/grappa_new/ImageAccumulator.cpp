#include "ImageAccumulator.h"

#include <chrono>
#include <ismrmrd/ismrmrd.h>
#include <boost/range/algorithm/copy.hpp>

#include "Reconstruction.h"

#include "Node.h"

#include "log.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;

    class AccumulationBuffer {
    public:
        AccumulationBuffer(Context context, OutputChannel &output);

        const Context context;
        OutputChannel &output;

        struct {
            std::vector<unsigned long> buffer_dimensions;
            int line_offset;
        };

        std::vector<hoNDArray<std::complex<float>>> image_buffers;

        void offer_acquisition(const Acquisition &acquisition);

        void copy_acquisition_to_buffer(const Acquisition &acquisition);
        void emit_reconstruction_job(const Acquisition &acquisition);

    private:
        hoNDArray<std::complex<float>> create_buffer();
        size_t get_receiver_channels(size_t default_receiver_channels = 1);
    };

    AccumulationBuffer::AccumulationBuffer(
            Context context,
            OutputChannel &output
    ) : context(std::move(context)), output(output) {

        ISMRMRD::EncodingSpace  r_space  = context.header.encoding[0].reconSpace;
        ISMRMRD::EncodingSpace  e_space  = context.header.encoding[0].encodedSpace;
        ISMRMRD::EncodingLimits e_limits = context.header.encoding[0].encodingLimits;

        line_offset = (e_space.matrixSize.y / 2) - e_limits.kspace_encoding_step_1->center;

        buffer_dimensions = {
            r_space.matrixSize.x,
            r_space.matrixSize.y,
            r_space.matrixSize.z,
            get_receiver_channels()
        };

        auto slices = e_limits.slice ? e_limits.slice->maximum + 1u : 1u;
        image_buffers = std::vector<hoNDArray<std::complex<float>>>(slices, create_buffer());
    }

    void AccumulationBuffer::offer_acquisition(const Acquisition &acquisition) {

        copy_acquisition_to_buffer(acquisition);

        auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);

        if (header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE)) emit_reconstruction_job(acquisition);
    }

    void AccumulationBuffer::copy_acquisition_to_buffer(
            const Acquisition &acquisition
    ) {
        auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        const auto &data = std::get<hoNDArray<std::complex<float>>>(acquisition);

        auto current_line      = header.idx.kspace_encode_step_1 + line_offset;
        auto current_slice     = header.idx.slice;
        auto current_partition = header.idx.kspace_encode_step_2;
        auto samples           = header.number_of_samples;

        auto &buffer = image_buffers[current_slice];

        // Copy the acquisition data to the buffer for each channel.
        for (size_t channel = 0; channel < header.active_channels; channel++) {

            auto destination = &buffer(0, size_t(current_line), current_partition, channel);
            auto source      = &data(0, channel);

            std::copy_n(source, samples, destination);
        }
    }

    void AccumulationBuffer::emit_reconstruction_job(const Acquisition &acquisition) {

        auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        auto slice = header.idx.slice;

        Grappa::Image image{};

        image.data = std::move(image_buffers[slice]);
        image.meta.slice = slice;

        boost::copy(header.position, image.meta.position.begin());
        boost::copy(header.read_dir, image.meta.read_dir.begin());
        boost::copy(header.phase_dir, image.meta.phase_dir.begin());
        boost::copy(header.slice_dir, image.meta.slice_dir.begin());
        boost::copy(header.patient_table_position, image.meta.table_pos.begin());

        output.push(std::move(image));

        image_buffers[slice] = create_buffer();
    }

    hoNDArray<std::complex<float>> AccumulationBuffer::create_buffer() {

        hoNDArray<std::complex<float>> buffer{buffer_dimensions};
        std::fill(buffer.begin(), buffer.end(), std::complex<float>(0.0f, 0.0f));

        return buffer;
    }

    size_t AccumulationBuffer::get_receiver_channels(size_t default_receiver_channels) {

        if (!context.header.acquisitionSystemInformation) return default_receiver_channels;
        if (!context.header.acquisitionSystemInformation->receiverChannels) return default_receiver_channels;

        return *context.header.acquisitionSystemInformation->receiverChannels;
    }
}

namespace Gadgetron::Grappa {

    ImageAccumulator::ImageAccumulator(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode<Acquisition>(props), context(context) {}

    void ImageAccumulator::process(TypedInputChannel<Acquisition> &in, OutputChannel &out) {

        AccumulationBuffer buffer{context, out};

        for (const auto &acquisition : in) buffer.offer_acquisition(acquisition);
    }

    GADGETRON_GADGET_EXPORT(ImageAccumulator);
}
