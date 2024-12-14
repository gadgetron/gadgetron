#include "FlowPhaseSubtractionGadget.h"
#include <queue>

#ifdef USE_OMP
#include <omp.h>
#endif

namespace Gadgetron {

void FlowPhaseSubtractionGadget::process(Core::InputChannel<mrd::Image<std::complex<float>>>& in,
                                         Core::OutputChannel& out) {

    const auto e_limits = this->header.encoding[0].encoding_limits;
    const auto sets = e_limits.set ? e_limits.set->maximum + 1 : 1;

    if (sets > 2)
        throw std::runtime_error("Phase subtraction only implemented for two sets");

    if (sets < 2) {
        std::move(in.begin(), in.end(), out.begin());
        return;
    }

    std::map<int, std::queue<mrd::Image<std::complex<float>>>> queues;

    for (auto image : in) {
        queues[image.head.set.value_or(0)].emplace(image);

        if (queues[0].empty() || queues[1].empty())
            continue;

        auto image1 = std::move(queues[0].front());
        auto image2 = std::move(queues[1].front());
        queues[0].pop();
        queues[1].pop();

        if (image1.head.image_index != image2.head.image_index)
            throw std::runtime_error("Mismatch in input indices detected");
        if (image1.data.size() != image2.data.size())
            throw std::runtime_error("Images must have same number of elements");

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (long i = 0; i < (long)image2.data.size(); i++) {
            std::complex<float> tmp =
                std::polar((std::abs(image1.data[i]) + std::abs(image2.data[i])) / 2.0f, std::arg(image2.data[i]) - std::arg(image1.data[i]));
            image2.data[i] = tmp;
        }

        image2.head.set = 0;

        out.push(std::move(image2));
    }
}

GADGETRON_GADGET_EXPORT(FlowPhaseSubtractionGadget)
} // namespace Gadgetron
