#include "NoiseSummaryGadget.h"

#include <boost/filesystem.hpp>
#include <limits>
#include <fstream>
#include "Dependency.h"
#include "NoiseAdjustGadget.h"

using namespace boost::filesystem;

namespace Gadgetron
{
    NoiseSummaryGadget::NoiseSummaryGadget()
    {
        processed_in_close_ = false;
    }

    NoiseSummaryGadget::~NoiseSummaryGadget()
    = default;

    int NoiseSummaryGadget::close(unsigned long flags)
    {
        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

        if ( !processed_in_close_ ) {
            processed_in_close_ = true;
            

            // list the content in the noise dependency folder
            // declear the attributes
            auto message = new Gadgetron::GadgetContainerMessage<DependencyQuery::Dependency>();
            auto& dependencies = message->getObjectPtr()->dependencies;

            auto legacy_id = this->noise_file.value();

            const auto legacy_prefix = std::string("GadgetronNoiseCovarianceMatrix_");

            auto is_prefix = [](auto&& potential_prefix, auto&& str){
                if (str.size() < potential_prefix.size()) return false;
                auto res = std::mismatch(potential_prefix.begin(),potential_prefix.end(),str.begin());
                if (res.first != potential_prefix.end()) return false;
                return true;
            };

            if (is_prefix(legacy_prefix,legacy_id)){
                legacy_id = legacy_id.substr(legacy_prefix.size());
            }

           auto  noise_covariance_list = this->context.storage.measurment.fetch<NoiseCovariance>(legacy_id,"noise_covariance");


            if ( noise_covariance_list.empty()) {
                dependencies.append("status", "failed");
            } else {

                const auto noise_covariance = noise_covariance_list[0];
                const auto& noise_covariance_matrix = noise_covariance.noise_covariance_matrix;

                size_t coils = noise_covariance_matrix.get_size(0);
                
                //Collect stats
                float mean_sigma = 0.0;
                float max_sigma = 0.0;
                float min_sigma = std::numeric_limits<float>::max();
                
                for (size_t c = 0; c < coils; c++) {
                    float sigma = std::sqrt(std::real(noise_covariance_matrix[c*coils+c]));
                    mean_sigma += sigma;
                    max_sigma = std::max(sigma,max_sigma);
                    min_sigma = std::min(sigma,min_sigma);
                }

                GDEBUG("Min Sigma: %f\n", min_sigma);
                GDEBUG("Max Sigma: %f\n", max_sigma);
                GDEBUG("Mean Sigma: %f\n", mean_sigma);

                dependencies.append("noise_dwell_time_us",noise_covariance.noise_dwell_time_us);
                dependencies.append("min_sigma",min_sigma);
                dependencies.append("max_sigma",max_sigma);
                dependencies.append("mean_sigma",mean_sigma);
                dependencies.append("channels", static_cast<long>(coils));
                dependencies.append("status", "success");
            }

            this->next()->putq(message);
            
        }
        
        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(NoiseSummaryGadget)

} // namespace Gadgetron
