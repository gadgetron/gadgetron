#include "AcquisitionBuffer.h"

#include "Context.h"
#include "Channel.h"
#include "Types.h"

#include "hoNDArray.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;

    size_t get_receiver_channels(const Context &context, size_t default_receiver_channels = 1) {

        if (!context.header.acquisitionSystemInformation) return default_receiver_channels;
        if (!context.header.acquisitionSystemInformation->receiverChannels) return default_receiver_channels;

        return *context.header.acquisitionSystemInformation->receiverChannels;
    }

    hoNDArray<std::complex<float>> create_buffer(const std::vector<unsigned long> &dimensions) {

        hoNDArray<std::complex<float>> buffer{dimensions};
        std::fill(buffer.begin(), buffer.end(), std::complex<float>(0.0f, 0.0f));

        return buffer;
    }
}

namespace Gadgetron::Grappa {

    AcquisitionBuffer::AcquisitionBuffer(Context ctx) : context(std::move(ctx)) {

        if (context.header.encoding.size() != 1) {
            throw std::runtime_error(
                    "This gadget only supports one encoding space; found " +
                    std::to_string(context.header.encoding.size())
            );
        }

        auto r_space  = context.header.encoding[0].reconSpace;
        auto e_space  = context.header.encoding[0].encodedSpace;
        auto e_limits = context.header.encoding[0].encodingLimits;

        auto slices = e_limits.slice ? e_limits.slice->maximum + 1u : 1u;

        internals.line_offset = (e_space.matrixSize.y / 2) - e_limits.kspace_encoding_step_1->center;
        internals.buffer_dimensions = {
                r_space.matrixSize.x,
                r_space.matrixSize.y,
                r_space.matrixSize.z,
                get_receiver_channels(context)
        };

        buffers = std::vector<hoNDArray<std::complex<float>>>(slices, create_buffer(internals.buffer_dimensions));
    }

    void AcquisitionBuffer::add(const Acquisition &acquisition) {

        auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        const auto &data = std::get<hoNDArray<std::complex<float>>>(acquisition);

        auto current_slice = header.idx.slice;
        auto current_line = header.idx.kspace_encode_step_1 + internals.line_offset;
        auto current_partition = header.idx.kspace_encode_step_2;
        auto samples = header.number_of_samples;

        auto &buffer = buffers[current_slice];

        // Copy the acquisition data to the buffer for each channel.
        for (size_t channel = 0; channel < header.active_channels; channel++) {

            auto destination = &buffer(0, size_t(current_line), current_partition, channel);
            auto source = &data(0, channel);

            std::copy_n(source, samples, destination);
        }

        for (const auto &hook : acquisition_hooks) hook(acquisition);
    }

    hoNDArray<std::complex<float>> AcquisitionBuffer::take(size_t index) {
        auto buffer = std::move(buffers[index]);
        clear(index);
        return buffer;
    }

    const hoNDArray<std::complex<float>> &AcquisitionBuffer::view(size_t index) {
        return buffers[index];
    }

    void AcquisitionBuffer::clear(size_t index) {
        buffers[index] = create_buffer(internals.buffer_dimensions);
    }

    void AcquisitionBuffer::add_acquisition_hook(std::function<void(const Core::Acquisition &)> fn) {
        acquisition_hooks.emplace_back(std::move(fn));
    }
}
