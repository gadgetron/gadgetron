#pragma once

#include "GadgetronTimer.h"
#include "Node.h"
#include "Types.h"
#include "gadgetron_mricore_export.h"
#include "hoNDArray.h"

#include <boost/filesystem/path.hpp>
#include <complex>
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

namespace Gadgetron {

        struct NoiseCovariance {
            ISMRMRD::IsmrmrdHeader header;
            float noise_dwell_time_us;
            hoNDArray<std::complex<float>> noise_covariance_matrix;
        };


        struct NoiseGatherer {
            hoNDArray<std::complex<float>> tmp_covariance;
            size_t number_of_samples = 0;
            float noise_dwell_time_us=0;
        };

        struct Prewhitener {
            hoNDArray<std::complex<float>> prewhitening_matrix;
        };

        struct LoadedNoise {
            hoNDArray<std::complex<float>> covariance;
            float noise_dwell_time_us;
        };

        struct IgnoringNoise {};
    class NoiseAdjustGadget : public Core::ChannelGadget<Core::Acquisition> {
    public:
        NoiseAdjustGadget(const Core::Context& context, const Core::GadgetProperties& props);

        void process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) override;


        using NoiseHandler = Core::variant<NoiseGatherer, LoadedNoise, Prewhitener, IgnoringNoise>;

    protected:
        NODE_PROPERTY(
            noise_dependency_prefix, std::string, "Prefix of noise dependency file", "GadgetronNoiseCovarianceMatrix");
        NODE_PROPERTY(perform_noise_adjust, bool, "Whether to actually perform the noise adjust", true);
        NODE_PROPERTY(pass_nonconformant_data, bool, "Whether to pass data that does not conform", true);
        NODE_PROPERTY(noise_dwell_time_us_preset, float, "Preset dwell time for noise measurement", 0.0);
        NODE_PROPERTY(
            scale_only_channels_by_name, std::string, "List of named channels that should only be scaled", "");
        NODE_PROPERTY(noise_dependency_folder, boost::filesystem::path, "Path to the working directory",
            boost::filesystem::temp_directory_path() / "gadgetron");

        const float receiver_noise_bandwidth;

        const std::string measurement_id;
        std::vector<size_t> scale_only_channels;

        // We will store/load a copy of the noise scans XML header to enable us to check which coil layout, etc.
        const ISMRMRD::IsmrmrdHeader current_ismrmrd_header;


        NoiseHandler noisehandler = IgnoringNoise{};





        template<class NOISEHANDLER>
        void add_noise(NOISEHANDLER& nh, const Core::Acquisition&) const ;

        template<class NOISEHANDLER>
        NoiseHandler handle_acquisition(NOISEHANDLER nh, Core::Acquisition&) const;




        template<class NOISEHANDLER>
        void save_noisedata(NOISEHANDLER& nh) const;

        NoiseHandler load_or_gather() const;
            
        bool first_run_;
    };
}
