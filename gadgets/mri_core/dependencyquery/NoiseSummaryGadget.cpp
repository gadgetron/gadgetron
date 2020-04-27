
#include "NoiseSummaryGadget.h"

#include "Dependency.h"
#include "NoiseAdjustGadget.h"
#include <boost/filesystem.hpp>
#include <fstream>
#include <limits>

using namespace boost::filesystem;

namespace Gadgetron {

    void NoiseSummaryGadget::process(Core::InputChannel<void>& input, Core::OutputChannel& output) {

        auto message       = DependencyQuery::Dependency{};
        auto& dependencies = message.dependencies;

        bool exists = context.storage.noise.contains(noise_file);
        dependencies.append("status", exists ? "success" : "failure");
        if (exists){

            auto noise_data = context.storage.noise.fetch<NoiseCovariance>(noise_file);

            size_t coils = noise_data.noise_covariance_matrix.get_size(0);

            // Collect stats
            float mean_sigma = 0.0;
            float max_sigma  = 0.0;
            float min_sigma  = std::numeric_limits<float>::max();

            for (size_t c = 0; c < coils; c++) {
                float sigma = std::sqrt(std::real(noise_data.noise_covariance_matrix[c * coils + c]));
                mean_sigma += sigma;
                if (sigma > max_sigma)
                    max_sigma = sigma;
                if (sigma < min_sigma)
                    min_sigma = sigma;
            }

            GDEBUG("Min Sigma: %f\n", min_sigma);
            GDEBUG("Max Sigma: %f\n", max_sigma);
            GDEBUG("Mean Sigma: %f\n", mean_sigma);

            dependencies.append("noise_dwell_time_us", noise_data.noise_dwell_time_us);
            dependencies.append("min_sigma", min_sigma);
            dependencies.append("max_sigma", max_sigma);
            dependencies.append("mean_sigma", mean_sigma);
            dependencies.append("channels", static_cast<long>(coils));
            dependencies.append("status", "success");
        }

        output.push(std::move(message));

    }
    NoiseSummaryGadget::NoiseSummaryGadget(const Core::Context& context, const Core::GadgetProperties& props)
        : Core::ChannelGadget<void>(context, props), context{ context } { }

} // namespace Gadgetron
