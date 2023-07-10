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
#include <nlohmann/json.hpp>
#include "sfndam.h"

using json = nlohmann::json;

namespace Gadgetron {
    
        struct NoiseCovariance {
            void SerializeToSfndam(std::ostream& stream) const
            {
                sfndam::sfndam<std::complex<float>> sf;
                json j;
                j["channels"] = channels_;
                j["sample_count"] = sample_count_;
                j["labels"] = labels_;
                j["noise_dwell_time_us"] = noise_dwell_time_us_;
                j["receiver_noise_bandwidth"] = receiver_noise_bandwidth_;
                sf.meta = j.dump();
                sf.array_dimensions = {static_cast<uint32_t>(channels_), static_cast<uint32_t>(channels_)};
                sf.data.resize(matrix_.get_number_of_elements());
                for (size_t i = 0; i < matrix_.get_number_of_elements(); i++) {
                    sf.data[i] = matrix_(i);
                }
                sfndam::serialize(sf, stream);
            }

            static NoiseCovariance DeserializeSfndam(std::istream& stream)
            {
                sfndam::sfndam<std::complex<float>> sf = sfndam::deserialize<std::complex<float>>(stream);
                json j = json::parse(sf.meta);
                auto labels = j["labels"].get<std::vector<std::string>>();
                NoiseCovariance out;
                out.labels_ = labels;
                out.channels_ = j["channels"].get<size_t>();
                out.sample_count_ = j["sample_count"].get<size_t>();
                out.noise_dwell_time_us_ = j["noise_dwell_time_us"].get<float>();
                out.receiver_noise_bandwidth_ = j["receiver_noise_bandwidth"].get<float>();
                out.matrix_ = hoNDArray<std::complex<float>(out.channels_, out.channels_);
                for (size_t i = 0; i < sf.data.size(); i++) {
                    out.matrix_(i) = sf.data[i];
                }
                return out;
            }

            size_t channels_;
            std::vector<std::string> labels_;
            hoNDArray<std::complex<float>> matrix_;
            size_t sample_count_;
            float noise_dwell_time_us_;
            float receiver_noise_bandwidth_;
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
        NoiseHandler handle_acquisition(NOISEHANDLER nh, Core::Acquisition&);


        Core::optional<NoiseCovariance> load_noisedata(const std::string& measurement_id) const;

        template<class NOISEHANDLER>
        void save_noisedata(NOISEHANDLER& nh);


        NoiseHandler load_or_gather() const;
        std::shared_ptr<MeasurementSpace> measurement_storage;

        // File name for file storage and retrival of noise covariance
        std::string noise_covariance_file_name = "";
    };
}
BOOST_HANA_ADAPT_STRUCT(Gadgetron::NoiseCovariance,header,noise_dwell_time_us,noise_covariance_matrix);
