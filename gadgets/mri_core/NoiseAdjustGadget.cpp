#include "NoiseAdjustGadget.h"
#include "hoArmadillo.h"
#include "hoMatrix.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_reductions.h"
#include "log.h"
#include "mri_core_utility.h"

#include <mrd/binary/protocols.h>

#include <boost/iterator/counting_iterator.hpp>
#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <typeinfo>

using namespace std::string_literals;
namespace bf = boost::filesystem;

namespace Gadgetron {
    namespace {

        float bandwidth_from_header(const mrd::Header& header) {
            if (header.acquisition_system_information && header.acquisition_system_information->relative_receiver_noise_bandwidth) {
                return header.acquisition_system_information->relative_receiver_noise_bandwidth.value();
            }
            return 0.793f;
        }

        std::string measurement_id_from_header(const mrd::Header& header) {
            if (header.measurement_information && header.measurement_information->measurement_id) {
                return header.measurement_information->measurement_id.value();
            }
            return ""s;
        }

        void normalize_covariance(NoiseGatherer& ng){
            if (ng.total_number_of_samples > 1) {
                ng.tmp_covariance /= std::complex<float>(ng.total_number_of_samples - 1);
                ng.normalized_number_of_samples = 1;
            }
        }

        // compare coil labels of noise and data
        // if number of channels are different, return false and order.size()==0
        // if any channels in noise cannot be found in data, return false and order.size()==0
        // if all channels in noise exist in data, but order is incorrect, return false, but  and order.size()==CHA
        // if all channels in nosie match channels in data, return true
        // order gives the matching order for src and dst coils
        // e.g. [2 1 0 3] means coil 0 of src matches coil 2 of dst etc.

        bool compare_coil_label(const std::vector<std::string>& src_coils,
            const std::vector<std::string>& dst_coils, std::vector<size_t>& order_in_src) {
            auto coil_name_comparer
                = [](const auto& coil1, const auto& coil2) { return coil1 == coil2; };
            bool labels_match = std::equal(
                src_coils.begin(), src_coils.end(), dst_coils.begin(), dst_coils.end(), coil_name_comparer);

            if (labels_match)
                return labels_match;

            if (!is_permutation(
                    src_coils.begin(), src_coils.end(), dst_coils.begin(), dst_coils.end(), coil_name_comparer))
                return false;

            order_in_src = std::vector<size_t>(dst_coils.size(), 0);
            std::iota(order_in_src.begin(), order_in_src.end(), 0);

            for (size_t d = 0u; d < order_in_src.size(); d++) {
                if (coil_name_comparer(dst_coils[d], src_coils[d]))
                    continue;
                auto coil_it    = std::find_if(src_coils.begin(), src_coils.end(),
                    [&](const auto& coil) { return coil_name_comparer(coil, dst_coils[d]); });
                order_in_src[d] = std::distance(src_coils.begin(), coil_it);
            }

            return labels_match;
        }

        hoNDArray<std::complex<float>> mask_channels(
            hoNDArray<std::complex<float>> noise_prewhitener_matrix, const std::vector<size_t>& scale_only_channels) {
            // Mask out scale  only channels
            size_t c                  = noise_prewhitener_matrix.get_size(0);
            std::complex<float>* dptr = noise_prewhitener_matrix.data();
            for (auto ch : scale_only_channels) {
                for (size_t i = 0; i < c; i++) {
                    for (size_t j = 0; j < c; j++) {
                        if ((i == ch || j == ch) && (i != j)) { // zero if scale only and not on diagonal
                            dptr[i * c + j] = std::complex<float>(0.0, 0.0);
                        }
                    }
                }
            }
            return std::move(noise_prewhitener_matrix);
        }

        hoNDArray<std::complex<float>> computeNoisePrewhitener(
            const hoNDArray<std::complex<float>>& noise_covariance_matrix) {

            auto noise_prewhitener_matrix = noise_covariance_matrix;
            size_t c                      = noise_prewhitener_matrix.get_size(0);
            float v                       = Gadgetron::asum(noise_covariance_matrix);
            if (v <= 0) {
                GDEBUG("Accumulated noise prewhitener is empty\n");
                for (size_t cha = 0; cha < c; cha++) {
                    noise_prewhitener_matrix(cha, cha) = 1;
                }
            } else {
                // Cholesky and invert lower triangular
                arma::cx_fmat noise_covf = as_arma_matrix(noise_prewhitener_matrix);
                noise_covf               = arma::inv(arma::trimatu(arma::chol(noise_covf)));
            }

            return noise_prewhitener_matrix;
        }

        std::vector<size_t> find_scale_only_channels(
            const std::string& scale_only_channels_by_name, const std::vector<mrd::CoilLabelType>& coillabels) {
            if (scale_only_channels_by_name.empty())
                return {};
            // Let's figure out if some channels are "scale_only"
            const std::string& uncomb_str = scale_only_channels_by_name;
            GDEBUG("SCALE ONLY: %s\n", uncomb_str.c_str());
            std::vector<std::string> uncomb;
            boost::split(uncomb, uncomb_str, boost::is_any_of(","));
            std::vector<size_t> scale_only_channels;

            for (unsigned int i = 0; i < uncomb.size(); i++) {
                std::string ch = boost::algorithm::trim_copy(uncomb[i]);
                if (std::find_if(
                        coillabels.begin(), coillabels.end(), [&](const auto& coil) { return ch == coil.coil_name; })
                    != coillabels.end())
                    scale_only_channels.push_back(i);
            }
            return scale_only_channels;
        }

        hoNDArray<std::complex<float>> reorder_noise_channels(
            hoNDArray<std::complex<float>> noise_covariance, const std::vector<size_t>& coil_order) {
            using namespace Indexing;
            // check whether to switch channel order
            auto CHA = noise_covariance.get_size(0);
            if ((coil_order.size() != CHA)
                || std::equal(coil_order.begin(), coil_order.end(), boost::counting_iterator<size_t>(0)))
                return std::move(noise_covariance);

            GDEBUG_STREAM("Require to reorder the noise covariance matrix to match the data ... ");
            hoNDArray<std::complex<float>> noise_covariance_reordered = noise_covariance;

            // switch row
            for (size_t n = 0; n < CHA; n++) {
                noise_covariance_reordered(n, slice) = noise_covariance(coil_order[n], slice);
            }

            // switch column
            for (size_t m = 0; m < CHA; m++) {
                noise_covariance(slice, m) = noise_covariance_reordered(slice, coil_order[m]);
            }

            return std::move(noise_covariance);
        }

        float calculate_scale_factor(
            float acquisition_dwell_time_us, float noise_dwell_time_us, float receiver_noise_bandwidth) {
            float noise_bw_scale_factor;
            if ((noise_dwell_time_us == 0.0f) || (acquisition_dwell_time_us == 0.0f)) {
                noise_bw_scale_factor = 1.0f;
            } else {
                noise_bw_scale_factor
                    = std::sqrt(2.0f * acquisition_dwell_time_us / noise_dwell_time_us * receiver_noise_bandwidth);
            }
            return noise_bw_scale_factor;
        }
    }

    NoiseAdjustGadget::NoiseAdjustGadget(const Core::Context& context, const Core::GadgetProperties& props)
        : Core::ChannelGadget<mrd::Acquisition>(context, props)
        , current_mrd_header(context.header)
        , receiver_noise_bandwidth{ bandwidth_from_header(context.header) }
        , measurement_id{ measurement_id_from_header(context.header) }
    {

        if (!perform_noise_adjust)
            return;

        GDEBUG("perform_noise_adjust_ is %d\n", perform_noise_adjust);
        GDEBUG("pass_nonconformant_data_ is %d\n", pass_nonconformant_data);
        GDEBUG("receiver_noise_bandwidth_ is %f\n", receiver_noise_bandwidth);

#ifdef USE_OMP
        omp_set_num_threads(1);
#endif // USE_OMP

        if (context.parameters.find("noisecovariancein") != context.parameters.end()) {
            noise_covariance_in = context.parameters.at("noisecovariancein");
            GDEBUG_STREAM("Input noise covariance matrix is provided as a parameter: " << noise_covariance_in);
        }

        if (context.parameters.find("noisecovarianceout") != context.parameters.end()) {
            noise_covariance_out = context.parameters.at("noisecovarianceout");
            GDEBUG_STREAM("Output noise covariance matrix is provided as a parameter: " << noise_covariance_out);
        }

        noisehandler = load_or_gather();
    }

    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::load_or_gather() const {
        auto noise_covariance = load_noisedata();

        if (noise_covariance) {
            size_t CHA = noise_covariance->matrix.get_size(0);
            if (noise_covariance->coil_labels.size() == CHA) {
                std::vector<std::string> current_coil_labels;
                if (current_mrd_header.acquisition_system_information) {
                    for (auto& l : current_mrd_header.acquisition_system_information->coil_label) {
                        current_coil_labels.push_back(l.coil_name);
                    }
                }

                std::vector<std::string> loaded_coil_labels;
                for (auto& l : noise_covariance->coil_labels) {
                    loaded_coil_labels.push_back(l.coil_name);
                }

                std::vector<size_t> coil_order_of_data_in_noise;
                bool labels_match = compare_coil_label(loaded_coil_labels,
                    current_coil_labels, coil_order_of_data_in_noise);

                if (!labels_match) {
                    // if number of channels in noise is different than data
                    // or
                    // if any channels in noise do not exist in data
                    if (CHA != current_coil_labels.size()) {
                        GDEBUG("Noise and measurement have different number of coils\n");
                    } else {
                        if (coil_order_of_data_in_noise.size() == CHA) {
                            GWARN_STREAM("Noise and measurement have different coils, but will be reordered ... ");
                            noise_covariance->matrix = reorder_noise_channels(noise_covariance->matrix, coil_order_of_data_in_noise);
                        } else {
                            GWARN_STREAM("Noise and measurement have different coils and cannot be reordered ... ");
                        }
                    }
                }
                return LoadedNoise{noise_covariance->matrix, noise_covariance->noise_dwell_time_us};

            } else if (current_mrd_header.acquisition_system_information) {
                GERROR("Noise covariance matrix is malformed. Number of labels does not match number of channels.");
            }
        }

        // No noise data found, gather it
        return NoiseGatherer{};
    }

    static bool is_noise(const mrd::Acquisition& acq) {
        return acq.head.flags.HasFlags(mrd::AcquisitionFlags::kIsNoiseMeasurement);
    }

    template <class NOISEHANDLER>
    void NoiseAdjustGadget::add_noise(NOISEHANDLER& nh, const mrd::Acquisition&) const {
    }

    template <> void NoiseAdjustGadget::add_noise(NoiseGatherer& ng, const mrd::Acquisition& acq) const {
        if (ng.tmp_covariance.empty()) {
            auto channels = acq.Coils();
            ng.tmp_covariance = hoNDArray<std::complex<float>>(channels, channels);
            std::fill(ng.tmp_covariance.begin(), ng.tmp_covariance.end(), std::complex<float>(0));
        }

        if (ng.noise_dwell_time_us == 0) {
            ng.noise_dwell_time_us = acq.head.sample_time_us.value_or(0);
        }

        auto dataM = as_arma_matrix(acq.data);
        auto covariance = as_arma_matrix(ng.tmp_covariance);
        covariance += dataM.t()*dataM;

        ng.total_number_of_samples += acq.Samples();
    }

    template <> void NoiseAdjustGadget::add_noise(NoiseHandler& nh, const mrd::Acquisition& acq) const {
        std::visit([&](auto& var) { this->add_noise(var, acq); }, nh);
    }
    template <class NOISEHANDLER> void NoiseAdjustGadget::save_noisedata(NOISEHANDLER& nh) {}

    template <> void NoiseAdjustGadget::save_noisedata(NoiseGatherer& ng) {
        if (ng.tmp_covariance.empty())
            return;

        normalize_covariance(ng);

        std::vector<mrd::CoilLabelType> coil_labels;
        for (auto& label : current_mrd_header.acquisition_system_information->coil_label) {
            coil_labels.push_back(label);
        }

        mrd::NoiseCovariance noise_covariance;
        noise_covariance.coil_labels = coil_labels;
        noise_covariance.sample_count = ng.total_number_of_samples;
        noise_covariance.noise_dwell_time_us = ng.noise_dwell_time_us;
        noise_covariance.receiver_noise_bandwidth = receiver_noise_bandwidth;
        noise_covariance.matrix = ng.tmp_covariance;

        if (!noise_covariance_out.empty()) {
            std::ofstream os(noise_covariance_out, std::ios::out | std::ios::binary);
            if (os.is_open()) {
                GDEBUG("Writing noise covariance to %s\n", noise_covariance_out.c_str());
                mrd::binary::MrdNoiseCovarianceWriter writer(os);
                writer.WriteNoiseCovariance(noise_covariance);
                writer.Close();
                os.flush();
                os.close();
            } else {
                GERROR("Unable to open file %s for writing noise covariance\n", noise_covariance_out.c_str());
                throw std::runtime_error("Unable to open file for writing noise covariance");
            }
        } else {
            GERROR_STREAM("Unable to save noise covariance. Noise covariance output file must be provided as a parameter");
            // throw std::runtime_error("Noise covariance output file must be provided as a parameter");
        }
    }

    template <> void NoiseAdjustGadget::save_noisedata(NoiseHandler& nh) {
        std::visit([&](auto& var) { this->save_noisedata(var); }, nh);
    }


    template <class NH>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(NH nh, mrd::Acquisition& acq) {
        return std::move(nh);
    };

    template <>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(
        Prewhitener pw, mrd::Acquisition& acq) {

        if (acq.Coils() == pw.prewhitening_matrix.get_size(0)) {
            auto dataM = as_arma_matrix(acq.data);
            auto pwm = as_arma_matrix(pw.prewhitening_matrix);
            dataM *= pwm;
        } else if (!this->pass_nonconformant_data) {
            throw std::runtime_error("Input data has different number of channels from noise data");
        }
        return std::move(pw);
    }

    template <>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(
        NoiseGatherer ng, mrd::Acquisition& acq) {
        if (ng.total_number_of_samples == 0) {
            return std::move(ng);
        }

        this->save_noisedata(ng);

        auto masked_covariance = mask_channels(ng.tmp_covariance, scale_only_channels);

        auto prewhitening_matrix = computeNoisePrewhitener(masked_covariance);
        prewhitening_matrix
            *= calculate_scale_factor(acq.head.sample_time_us.value_or(0), ng.noise_dwell_time_us, receiver_noise_bandwidth);
        return handle_acquisition(Prewhitener{ prewhitening_matrix }, acq);
    }

    template <>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(
        LoadedNoise ln, mrd::Acquisition& acq)  {
        auto masked_covariance   = mask_channels(std::move(ln.covariance), scale_only_channels);
        auto prewhitening_matrix = computeNoisePrewhitener(masked_covariance);
        prewhitening_matrix
            *= calculate_scale_factor(acq.head.sample_time_us.value_or(0), ln.noise_dwell_time_us, receiver_noise_bandwidth);
        return handle_acquisition(Prewhitener{ prewhitening_matrix }, acq);
    }

    template <>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(
        NoiseHandler nh, mrd::Acquisition& acq) {
        return std::visit([&](auto var) { return this->handle_acquisition<decltype(var)>(std::move(var), acq); }, std::move(nh));
    }

    void NoiseAdjustGadget::process(Core::InputChannel<mrd::Acquisition>& input, Core::OutputChannel& output) {

        scale_only_channels = current_mrd_header.acquisition_system_information
                                  ? find_scale_only_channels(scale_only_channels_by_name,
                                      current_mrd_header.acquisition_system_information->coil_label)
                                  : std::vector<size_t>{};

        for (auto acq : input) {
            if (is_noise(acq)) {
                add_noise(noisehandler, acq);
                continue;
            }
            noisehandler = handle_acquisition(std::move(noisehandler), acq);
            output.push(std::move(acq));
        }

        this->save_noisedata(noisehandler);
    }

    /** Returns NoiseCovariance if loaded from file/stream, otherwise None */
    std::optional<mrd::NoiseCovariance> NoiseAdjustGadget::load_noisedata() const {
        if (!noise_covariance_in.empty()) {
            std::ifstream file(noise_covariance_in, std::ios::binary);
            if (!file) {
                GERROR("Could not open noise covariance file %s\n", noise_covariance_in.c_str());
                GWARN("Falling back to noise gathering\n");
                // throw std::runtime_error("Could not open noise covariance file");
                return std::nullopt;
            }
            mrd::binary::MrdNoiseCovarianceReader reader(file);
            mrd::NoiseCovariance noise_covariance;
            reader.ReadNoiseCovariance(noise_covariance);
            reader.Close();
            file.close();
            return noise_covariance;
        }

        return std::nullopt;
    }

    GADGETRON_GADGET_EXPORT(NoiseAdjustGadget)

} // namespace Gadgetron
