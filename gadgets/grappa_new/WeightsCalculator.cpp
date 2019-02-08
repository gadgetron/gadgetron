#include "WeightsCalculator.h"

#include <set>

#include "common/AcquisitionBuffer.h"
#include "common/slice.h"

#include "SliceAccumulator.h"
#include "Unmixing.h"

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

        void operator()(const Acquisition &acquisition) {
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

        void operator()(const Acquisition &acquisition) {

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

        void operator()(const Acquisition &acquisition) {

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

    class WeightCore {
    public:
        WeightCore() {

        }

        hoNDArray<std::complex<float>> calculate_weights(
                hoNDArray<std::complex<float>> data,
                std::array<size_t, 4> region_of_support,
                size_t acceleration_factor
        );

        hoNDArray<std::complex<float>> calculate_coil_map(

        );

    private:
        // In an attempt to avoid unnecessary allocations, we maintain (and overwrite) some buffers
        // (coil map and image kernels) between weight calculations.

        size_t target_channels; // Initialize from context? Yes, everything in context/props.
    };

    hoNDArray<std::complex<float>> WeightCore::calculate_coil_map(

    ) {
        return hoNDArray<std::complex<float>>();
    }

    hoNDArray<std::complex<float>> WeightCore::calculate_weights(
            hoNDArray<std::complex<float>> data,
            std::array<size_t, 4> region_of_support,
            size_t acceleration_factor
    ) {
        GINFO_STREAM("Acceleration Factor: " << acceleration_factor);
        GINFO_STREAM("Region of support:");
        GINFO_STREAM("\tStart R0: " << region_of_support[0]);
        GINFO_STREAM("\t  End R0: " << region_of_support[1]);
        GINFO_STREAM("\tStart E1: " << region_of_support[2]);
        GINFO_STREAM("\t  End E1: " << region_of_support[3]);

        return hoNDArray<std::complex<float>>();

        // Reorder the Channels - Combined channels, then uncombined channels after.



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

        std::set<size_t> updated_slices{};

        AcquisitionBuffer buffer{context};
        SupportMonitor support_monitor{context};
        AccelerationMonitor acceleration_monitor{context};

        buffer.add_post_update_callback([&](auto &acq) { updated_slices.insert(slice_of(acq)); });
        buffer.add_post_update_callback([&](auto &acq) { acceleration_monitor(acq); });
        buffer.add_post_update_callback([&](auto &acq) { support_monitor(acq); });
        buffer.add_pre_update_callback(DirectionMonitor{buffer, support_monitor, acceleration_monitor});

        while (true) {
            auto slices = take_available_slices(in);
            buffer.add(slices);

            for (auto index : updated_slices) {
//                Grappa::Weights weights {
//                        { index },
//                        calculate_weights(
//                                buffer.view(index),
//                                support_monitor.region_of_support(index),
//                                acceleration_monitor.acceleration_factor(index)
//                        )
//                };
//                out.push(std::move(weights));
            }

            updated_slices.clear();
        }
    }

    GADGETRON_GADGET_EXPORT(WeightsCalculator);
}