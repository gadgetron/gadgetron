#include "WeightsCalculator.h"

#include "common/AcquisitionBuffer.h"
#include "SliceAccumulator.h"
#include "Reconstruction.h"

#include "Gadget.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;

    std::vector<Slice> take_available_slices(TypedInputChannel<Slice> &input) {

        std::vector<Slice> slices{};

        slices.emplace_back(input.pop());

        while(auto opt_slice = input.try_pop()) {
            slices.emplace_back(std::move(*opt_slice));
        }

        return std::move(slices);
    }

    class DirectionMonitor {
    public:
        explicit DirectionMonitor(AcquisitionBuffer &buffer) : buffer(buffer) {}

        void operator()(const Acquisition &acquisition) {

        }
    private:
        AcquisitionBuffer &buffer;

        std::array<float, 3> position, read_dir, phase_dir, slice_dir;
    };

    class AccelerationMonitor {
    public:
        explicit AccelerationMonitor(AcquisitionBuffer &buffer) : buffer(buffer) {}

        void operator()(const Acquisition &acquisition) {

        }
    private:
        AcquisitionBuffer &buffer;
    };
}

namespace Gadgetron::Grappa {

    WeightsCalculator::WeightsCalculator(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : TypedGadgetNode<Slice>(props), context(context) {}

    void WeightsCalculator::process(TypedInputChannel<Slice> &in, OutputChannel &out) {

        AcquisitionBuffer buffer{context};
        DirectionMonitor discard_buffer_on_scan_direction_change{buffer};
        AccelerationMonitor discard_buffer_on_acceleration_factory_change{buffer};

        buffer.add_acquisition_hook(discard_buffer_on_scan_direction_change);
        buffer.add_acquisition_hook(discard_buffer_on_acceleration_factory_change);

        while (true) {
            auto slices = take_available_slices(in);
            buffer.add(slices);

            // TODO: Recalculate the weights. Emit the weights.
        }
    }

    GADGETRON_GADGET_EXPORT(WeightsCalculator);
}