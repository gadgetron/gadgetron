#include "AcquisitionBuffer.h"

#include <map>
#include <set>

#include "Context.h"
#include "Channel.h"
#include "Types.h"

#include "hoNDArray.h"

#include "grappa_common.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;

    template<class T>
    std::set<T> all_values_in_range(T low, T high) {
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

        if (r_space.matrixSize.z != 1) {
            throw std::runtime_error("RT Grappa works only with 2D images. 3D output requested.");
        }

        internals.number_of_samples = r_space.matrixSize.x;
        internals.number_of_lines = r_space.matrixSize.y;

        internals.line_offset = (e_space.matrixSize.y / 2) - e_limits.kspace_encoding_step_1->center;
        internals.expected_lines = all_values_in_range<uint32_t>(
                e_limits.kspace_encoding_step_1->minimum,
                e_limits.kspace_encoding_step_1->maximum
        );

        buffers = std::map<size_t, buffer>{};
    }

    void AcquisitionBuffer::add(const AnnotatedAcquisition &acquisition) {

        for (auto &fn : pre_update_callbacks) fn(acquisition);

        auto &header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        const auto &data = std::get<hoNDArray<std::complex<float>>>(acquisition);

        if (header.idx.kspace_encode_step_2 != 0) {
            throw std::runtime_error("RT Grappa works only with 2D data. 3D data received.");
        }

        auto current_slice = header.idx.slice;
        auto current_line = header.idx.kspace_encode_step_1 + internals.line_offset;

        if (!buffers.count(current_slice)) {
            buffers[current_slice] = create_buffer({
                internals.number_of_samples,
                internals.number_of_lines,
                header.active_channels
            });
        }

        auto &buffer = buffers[current_slice];
        buffer.sampled_lines.insert(header.idx.kspace_encode_step_1);

        // Copy the acquisition data to the buffer for each channel.
        for (size_t channel = 0; channel < header.active_channels; channel++) {

            auto destination = &buffer.data(0, size_t(current_line), channel);
            auto source = &data(0, channel);

            std::copy_n(source, header.number_of_samples, destination);
        }

        for (auto &fn : post_update_callbacks) fn(acquisition);
    }

    hoNDArray<std::complex<float>> AcquisitionBuffer::take(size_t index) {
        auto buffer = buffers[index];
        clear(index);
        return std::move(buffer.data);
    }

    const hoNDArray<std::complex<float>> &AcquisitionBuffer::view(size_t index) const {
        return buffers.at(index).data;
    }

    void AcquisitionBuffer::clear(size_t index) {
        buffers.erase(index);
    }

    bool AcquisitionBuffer::is_fully_sampled(size_t index) const {
        return internals.expected_lines == buffers.at(index).sampled_lines;
    }

    AcquisitionBuffer::buffer AcquisitionBuffer::create_buffer(const std::vector<size_t> &dimensions) {

        buffer buffer {
            hoNDArray<std::complex<float>>(dimensions),
            {}
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
