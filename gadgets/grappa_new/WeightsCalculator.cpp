#include "WeightsCalculator.h"

#include <set>

#include "common/AcquisitionBuffer.h"
#include "common/grappa_common.h"

#include "SliceAccumulator.h"
#include "Unmixing.h"

#include "GadgetronTimer.h"
#include "hoNDFFT.h"
#include "hoNDArray_iterators.h"
#include "hoNDArray_utils.h"
#include "mri_core_grappa.h"
#include "mri_core_coil_map_estimation.h"

#include "Gadget.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;

    // A similar function should be available in the std library at some point.
    template <class T, std::size_t N>
    std::array<std::remove_cv_t<T>, N> to_array(T (&a)[N])
    {
        std::array<std::remove_cv_t<T>, N> array{};
        std::copy(std::begin(a), std::end(a), array.begin());
        return array;
    }

    std::vector<Slice> take_available_slices(TypedInputChannel<Slice> &input) {

        std::vector<Slice> slices{};

        slices.emplace_back(input.pop());

        while(auto opt_slice = input.try_pop()) {
            slices.emplace_back(std::move(*opt_slice));
        }

        GDEBUG_STREAM("Read " << slices.size() << " available slice(s).");

        return std::move(slices);
    }

    class SupportMonitor {
    public:
        explicit SupportMonitor(const Context &context) {
            auto e_limits = context.header.encoding[0].encodingLimits;
            auto slices = e_limits.slice ? e_limits.slice->maximum + 1u : 1u;

            regions = std::vector<std::array<size_t, 4>>(slices, clear_region);
        }

        void operator()(const AnnotatedAcquisition &acquisition) {
            auto old_region = regions[slice_of(acquisition)];
            regions[slice_of(acquisition)] = std::array<size_t, 4> {
                0,
                samples_in(acquisition) - 1,
                std::min(old_region[2], line_of(acquisition)),
                std::max(old_region[3], line_of(acquisition))
            };
        }

        std::array<size_t, 4> region_of_support(size_t slice) {
            return regions[slice];
        }

        void clear(size_t slice) {
            regions[slice] = clear_region;
        }

    private:
        std::vector<std::array<size_t, 4>> regions;

        static constexpr std::array<size_t, 4> clear_region = {
            0, 0, std::numeric_limits<size_t>::max(), 0
        };
    };

    class AccelerationMonitor {
    public:
        AccelerationMonitor(const Context &context) {
            auto e_limits = context.header.encoding[0].encodingLimits;
            auto slices = e_limits.slice ? e_limits.slice->maximum + 1u : 1u;

            previous_line = std::vector<optional<size_t>>(slices, none);
            acceleration = std::vector<optional<size_t>>(slices, none);
        }

        void operator()(const AnnotatedAcquisition &acquisition) {

            if(previous_line[slice_of(acquisition)]) {
                if (line_of(acquisition) < previous_line[slice_of(acquisition)].get()) {
                    acceleration[slice_of(acquisition)] = none;
                }
                else {
                    acceleration[slice_of(acquisition)] = line_of(acquisition) - previous_line[slice_of(acquisition)].get();
                }
            }

            previous_line[slice_of(acquisition)] = line_of(acquisition);
        }

        size_t acceleration_factor(size_t slice) {
            return acceleration[slice].get();
        }

        void clear(size_t slice) {
            previous_line[slice] = acceleration[slice] = none;
        }

    private:
        std::vector<optional<size_t>> previous_line;
        std::vector<optional<size_t>> acceleration;
    };

    class DirectionMonitor {
    public:
        explicit DirectionMonitor(AcquisitionBuffer &buffer, SupportMonitor &support, AccelerationMonitor &acceleration)
        : buffer(buffer), support(support), acceleration(acceleration) {
            position = read_dir = phase_dir = slice_dir = {0.0, 0.0, 0.0};
        }

        void operator()(const AnnotatedAcquisition &acquisition) {

            auto header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);

            if (position == to_array(header.position) &&
                read_dir == to_array(header.read_dir) &&
                phase_dir == to_array(header.phase_dir) &&
                slice_dir == to_array(header.slice_dir)) {
                return;
            }

            position = to_array(header.position);
            read_dir = to_array(header.read_dir);
            phase_dir = to_array(header.phase_dir);
            slice_dir = to_array(header.slice_dir);

            clear(slice_of(acquisition));
        }

        void clear(size_t slice) {
            buffer.clear(slice);
            support.clear(slice);
            acceleration.clear(slice);
        }


    private:
        AcquisitionBuffer &buffer;
        SupportMonitor &support;
        AccelerationMonitor &acceleration;

        std::array<float, 3> position, read_dir, phase_dir, slice_dir;
    };

    // ---------------------------------------------------------------------------- //
    // ---------------------------------------------------------------------------- //
    // ---------------------------------------------------------------------------- //

    class WeightsCore {
    public:
        hoNDArray<std::complex<float>> calculate_weights(
                hoNDArray<std::complex<float>> data,
                std::array<size_t, 4> region_of_support,
                size_t acceleration_factor,
                uint64_t n_combined_channels
        );

        const hoNDArray<std::complex<float>> &
        estimate_coil_map(
                const hoNDArray<std::complex<float>> &data
        );

        const hoNDArray<std::complex<float>>
        keep_combined_channels(
                uint64_t n_combined_channels,
                const hoNDArray<std::complex<float>> &data
        );

        struct {
            size_t ks, power;
        } coil_map_params;

        struct {
            size_t width, height;
            float threshold;
        } kernel_params;

        struct {
            // We maintain a few buffers to avoid reallocating them repeatedly.
            hoNDArray<std::complex<float>> image, coil_map, convolution_kernel, image_domain_kernel;
            hoNDArray<float> g_factor;
        } buffers;
    };

    const hoNDArray<std::complex<float>> &WeightsCore::estimate_coil_map(const hoNDArray<std::complex<float>> &data) {

        hoNDFFT<float>::instance()->ifft2c(data, buffers.image);
        Gadgetron::coil_map_2d_Inati(buffers.image, buffers.coil_map, coil_map_params.ks, coil_map_params.power);

        return buffers.coil_map;
    }

    const hoNDArray<std::complex<float>> WeightsCore::keep_combined_channels(
            uint64_t n_combined_channels,
            const hoNDArray<std::complex<float>> &data
    ) {
        return hoNDArray<std::complex<float>>(
            data.dimensions()[0],
            data.dimensions()[1],
            n_combined_channels,
            const_cast<std::complex<float> *>(data.data())
        );
    }

    hoNDArray<std::complex<float>> WeightsCore::calculate_weights(
            hoNDArray<std::complex<float>> data,
            std::array<size_t, 4> region_of_support,
            size_t acceleration_factor,
            uint64_t n_combined_channels
    ) {
        // TODO: Optimize accel_factor == 1;

        GadgetronTimer weights_timer("Grappa weights calculation");
        data.squeeze();

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t CHA = data.get_size(2);

        auto combined_data = keep_combined_channels(n_combined_channels, data);
        auto coil_map = estimate_coil_map(combined_data);

        Gadgetron::grappa2d_calib_convolution_kernel(
                data,
                combined_data,
                acceleration_factor,
                kernel_params.threshold,
                kernel_params.width,
                kernel_params.height,
                region_of_support[0],
                region_of_support[1],
                region_of_support[2],
                region_of_support[3],
                buffers.convolution_kernel
        );

        Gadgetron::grappa2d_image_domain_kernel(
                buffers.convolution_kernel,
                RO,
                E1,
                buffers.image_domain_kernel
        );

        hoNDArray<std::complex<float>> unmixing_coefficients(RO, E1, CHA);
        Gadgetron::grappa2d_unmixing_coeff(
                buffers.image_domain_kernel,
                coil_map,
                acceleration_factor,
                unmixing_coefficients,
                buffers.g_factor
        );

        return unmixing_coefficients;
    }

    // ---------------------------------------------------------------------------- //
    // ---------------------------------------------------------------------------- //
    // ---------------------------------------------------------------------------- //
}

namespace Gadgetron::Grappa {

    WeightsCalculator::WeightsCalculator(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode<Slice>(props), context(context) {}

    void WeightsCalculator::process(TypedInputChannel<Slice> &in, OutputChannel &out) {

        std::set<uint16_t> updated_slices{};
        uint16_t n_combined_channels = 0, n_uncombined_channels = 0;

        AcquisitionBuffer buffer{context};
        SupportMonitor support_monitor{context};
        AccelerationMonitor acceleration_monitor{context};

        buffer.add_post_update_callback([&](auto &acq) { updated_slices.insert(slice_of(acq)); });
        buffer.add_post_update_callback([&](auto &acq) { acceleration_monitor(acq); });
        buffer.add_post_update_callback([&](auto &acq) { support_monitor(acq); });
        buffer.add_post_update_callback([&](auto &acq) {
            n_combined_channels = combined_channels(acq);
            n_uncombined_channels = uncombined_channels(acq);
        });

        buffer.add_pre_update_callback(DirectionMonitor{buffer, support_monitor, acceleration_monitor});

        WeightsCore core{
                {coil_map_estimation_ks, coil_map_estimation_power},
                {block_size_samples, block_size_lines, convolution_kernel_threshold}
        };

        while (true) {
            auto slices = take_available_slices(in);
            buffer.add(slices);

            for (auto index : updated_slices) {
                Grappa::Weights weights {
                        { index, n_combined_channels, n_uncombined_channels },
                        core.calculate_weights(
                                buffer.view(index),
                                support_monitor.region_of_support(index),
                                acceleration_monitor.acceleration_factor(index),
                                n_combined_channels
                        )
                };
                out.push(std::move(weights));
            }

            updated_slices.clear();
        }
    }

    GADGETRON_GADGET_EXPORT(WeightsCalculator);
}