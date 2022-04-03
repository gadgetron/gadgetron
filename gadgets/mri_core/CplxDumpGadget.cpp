#include "CplxDumpGadget.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"

namespace Gadgetron {
void CplxDumpGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) {

    std::vector<Core::Acquisition> buffer;

    for (auto [header, acq, traj] : in) {
        if (ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(header.flags))
            continue;
        buffer.emplace_back(header, acq, traj);
        out.push(Core::Acquisition{header, std::move(acq), std::move(traj)});
    }

    GDEBUG("CplxDumpGadget::close...\n");
    GDEBUG("Number of items on Q: %d\n", buffer.size());

    if (buffer.empty())
        return;

    const auto& [header, data, traj] = buffer.back();
    auto dims = data.dimensions();
    auto dims_profile = data.dimensions();
    dims.push_back(buffer.size());

    // Allocate array for result
    hoNDArray<std::complex<float>> result(dims);

    for (size_t i = 0; i < buffer.size(); i++) {
        using namespace Gadgetron::Indexing;
        result(slice, slice, i) = std::get<1>(buffer[i]);
    }

    // Reshape to get the coil dimension as the last
    std::vector<size_t> order = {0, 2, 1};
    result = permute(result, order);

    // Write out the result
    write_nd_array<std::complex<float>>(&result, filename.c_str());
}

GADGETRON_GADGET_EXPORT(CplxDumpGadget)

} // namespace Gadgetron
