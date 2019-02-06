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

        GDEBUG_STREAM("Read " << slices.size() << " available slices.");

        return std::move(slices);
    }

    class SupportMonitor {
    public:
        explicit SupportMonitor(const Context &context) {
            auto e_limits = context.header.encoding[0].encodingLimits;
            auto slices = e_limits.slice ? e_limits.slice->maximum + 1u : 1u;

            regions = std::vector<std::pair<size_t, size_t>>(
                    slices,
                    std::make_pair(std::numeric_limits<size_t>::max(), size_t(0))
            );
        }

        void operator()(const Acquisition &acquisition) {
            auto old_region = regions[slice_of(acquisition)];
            regions[slice_of(acquisition)] = std::make_pair(
                    std::min(old_region.first, line_of(acquisition)),
                    std::max(old_region.second, line_of(acquisition))
            );
        }

        std::pair<size_t, size_t> region_of_support(size_t slice) {
            return regions[slice];
        }

        void clear(size_t slice) {
            regions[slice] = std::make_pair(
                    std::numeric_limits<size_t>::max(),
                    size_t(0)
            );
        }

    private:
        std::vector<std::pair<size_t, size_t>> regions;
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

    hoNDArray<std::complex<float>> calculate_weights(
            hoNDArray<std::complex<float>> data,
            std::pair<size_t, size_t> region_of_support,
            size_t acceleration_factor
    ) {
        GINFO_STREAM("Hello, I'm calculating some weights.");

        return hoNDArray<std::complex<float>>();
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
                Grappa::Weights weights {
                        { index },
                        calculate_weights(
                                buffer.view(index),
                                support_monitor.region_of_support(index),
                                acceleration_monitor.acceleration_factor(index)
                        )
                };
                out.push(std::move(weights));
            }

            updated_slices.clear();
        }
    }

    GADGETRON_GADGET_EXPORT(WeightsCalculator);
}