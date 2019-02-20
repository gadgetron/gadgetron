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

    template<class T>
    std::set<T> create_set(T low, T high) {
        std::set<T> s{};
        for (auto i = low; i <= high; i++) s.insert(i);
        return s;
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

        if (r_space.matrixSize.z != 1) {
            throw std::runtime_error("RT Grappa works only with 2D images. 3D output requested.");
        }

        internals.line_offset = (e_space.matrixSize.y / 2) - e_limits.kspace_encoding_step_1->center;
        internals.expected_lines = create_set(
                e_limits.kspace_encoding_step_1->minimum,
                e_limits.kspace_encoding_step_1->maximum
        );
        internals.buffer_dimensions = {
                r_space.matrixSize.x,
                r_space.matrixSize.y,
                get_receiver_channels(context)
        };

        buffers = std::vector<buffer>(slices, create_buffer());
    }

    void AcquisitionBuffer::add(const AnnotatedAcquisition &acquisition) {

        for (auto &fn : pre_update_callbacks) fn(acquisition);

        auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        const auto &data = std::get<hoNDArray<std::complex<float>>>(acquisition);

        if (header.idx.kspace_encode_step_2 != 0) {
            throw std::runtime_error("RT Grappa works only on 2D data. 3D data received.");
        }

        auto current_slice = header.idx.slice;
        auto current_line = header.idx.kspace_encode_step_1 + internals.line_offset;
        auto samples = header.number_of_samples;

        auto &buffer = buffers[current_slice];
        buffer.sampled_lines.insert(header.idx.kspace_encode_step_1);

        // Copy the acquisition data to the buffer for each channel.
        for (size_t channel = 0; channel < header.active_channels; channel++) {

            auto destination = &buffer.data(0, size_t(current_line), channel);
            auto source = &data(0, channel);

            std::copy_n(source, samples, destination);
        }

        for (auto &fn : post_update_callbacks) fn(acquisition);
    }

    hoNDArray<std::complex<float>> AcquisitionBuffer::take(size_t index) {
        auto buffer = std::move(buffers[index]);
        clear(index);
        return buffer.data;
    }

    const hoNDArray<std::complex<float>> &AcquisitionBuffer::view(size_t index) {
        return buffers[index].data;
    }

    void AcquisitionBuffer::clear(size_t index) {
        buffers[index] = create_buffer();
    }

    bool AcquisitionBuffer::is_fully_sampled(size_t index) {
        return internals.expected_lines == buffers[index].sampled_lines;
    }

    AcquisitionBuffer::buffer AcquisitionBuffer::create_buffer() {

        buffer buffer {
            hoNDArray<std::complex<float>>(internals.buffer_dimensions),
            std::set<uint16_t>()
        };

        std::fill(buffer.data.begin(), buffer.data.end(), std::complex<float>(0.0f, 0.0f));

        return buffer;
    }


    void AcquisitionBuffer::add_pre_update_callback(std::function<void(const AnnotatedAcquisition &)> fn) {
        pre_update_callbacks.emplace_back(std::move(fn));
    }

    void AcquisitionBuffer::add_post_update_callback(std::function<void(const AnnotatedAcquisition &)> fn) {
        post_update_callbacks.emplace_back(std::move(fn));
    }
}
