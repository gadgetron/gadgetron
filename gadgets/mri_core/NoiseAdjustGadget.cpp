#include "NoiseAdjustGadget.h"
#include "hoArmadillo.h"
#include "hoMatrix.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_reductions.h"
#include "io/primitives.h"
#include "io/ismrmrd_types.h"
#include "log.h"
#include <boost/iterator/counting_iterator.hpp>
#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cmath>
#include <iomanip>
#include <limits>
#include <typeinfo>

using namespace std::string_literals;
namespace bf = boost::filesystem;

namespace Gadgetron {
    namespace {

        template <class T> T value_or(const ISMRMRD::Optional<T>& opt, T default_value) {
            return opt ? *opt : default_value;
        }

        float bandwidth_from_header(const ISMRMRD::IsmrmrdHeader& header) {
            return value_or(header.acquisitionSystemInformation->relativeReceiverNoiseBandwidth, 0.793f);
        }

        void normalize_covariance(NoiseGatherer& ng){
            if (ng.total_number_of_samples > 1) {
                ng.tmp_covariance /= std::complex<float>(ng.total_number_of_samples - 1);
                ng.normalized_number_of_samples = 1;
            }
        }

        std::string to_string(const std::vector<ISMRMRD::CoilLabel>& coils) {
            std::stringstream sstream;
            for (auto i = 0u; i < coils.size(); i++)
                sstream << "Coil " << i << " - " << coils[i].coilNumber << " - " << coils[i].coilName << std::endl;
            return sstream.str();
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
            const std::string& scale_only_channels_by_name, const std::vector<ISMRMRD::CoilLabel>& coillabels) {
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
                        coillabels.begin(), coillabels.end(), [&](const auto& coil) { return ch == coil.coilName; })
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
                //hoNDArrayView<std::complex<float>,1,false> f = noise_covariance_reordered(n,slice);
                noise_covariance_reordered(n, slice) = noise_covariance(coil_order[n], slice);
            }

            // switch column
            for (size_t m = 0; m < CHA; m++) {
                noise_covariance(slice, m) = noise_covariance_reordered(slice, coil_order[m]);
            }

            return std::move(noise_covariance);
        }

        void print_covariance_matrix(const hoNDArray<std::complex<float>>& covariance, const std::vector<std::string>& labels = {}, size_t num_samples = 0) {
            const char* shades = " \xE2\x96\x91\xE2\x96\x92\xE2\x96\x93\xE2\x96\x88";
            const int shade_widths[] = {1, 3, 3, 3, 3}; // byte widths: ' '=1, UTF-8 block chars=3
            size_t n = covariance.get_size(0);

            // Compute absolute values and find maximum
            float max_val = 0.0f;
            for (size_t i = 0; i < n * n; i++) {
                float a = std::abs(covariance.data()[i]);
                if (a > max_val) max_val = a;
            }

            if (max_val == 0.0f) max_val = 1.0f;

            std::string output;
            output += "\033[32mNoise covariance matrix:\033[0m\n";
            std::ostringstream max_val_stream;
            max_val_stream << std::scientific << std::setprecision(3) << max_val;
            output += "Noise covariance matrix (" + std::to_string(n) + " x " + std::to_string(n) + " channels), max = " + max_val_stream.str() + (num_samples ? ", samples = " + std::to_string(num_samples) : "") + ":\n";
            size_t max_label_width = 0;
            for (const auto& label : labels) {
                max_label_width = std::max(max_label_width, label.size());
            }

            float min_nonzero_diag = std::numeric_limits<float>::infinity();
            for (size_t i = 0; i < n; i++) {
                float diag_abs = std::abs(covariance(i, i));
                if (diag_abs > 0.0f) {
                    min_nonzero_diag = std::min(min_nonzero_diag, diag_abs);
                }
            }

            int scale_power = 0;
            if (std::isfinite(min_nonzero_diag) && min_nonzero_diag < 1.0f) {
                scale_power = static_cast<int>(std::ceil(-std::log10(min_nonzero_diag)));
            }
            float scale_factor = std::pow(10.0f, static_cast<float>(scale_power));
            output += "diag scale factor = 10^" + std::to_string(scale_power) + "\n";

            for (size_t i = 0; i < n; i++) {
                if (i < labels.size()) {
                    output += labels[i];
                    output.append(max_label_width - labels[i].size(), ' ');
                    output += " | ";
                }
                for (size_t j = 0; j < n; j++) {
                    float normalized = std::abs(covariance(i, j)) / max_val;
                    int idx = std::min((int)(normalized * 5.0f), 4);
                    const char* p = shades;
                    for (int k = 0; k < idx; k++) p += shade_widths[k];
                    output.append(p, shade_widths[idx]);
                }
                std::ostringstream diagonal_stream;
                diagonal_stream << std::fixed << std::setprecision(2) << std::setw(5)
                                << std::abs(covariance(i, i)) * scale_factor;
                output += " | diag = " + diagonal_stream.str();
                output += '\n';
            }
            GDEBUG("%s", output.c_str());
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
        : Core::ChannelGadget<Core::Acquisition>(context, props)
        , current_ismrmrd_header(context.header)
        , receiver_noise_bandwidth{ bandwidth_from_header(context.header) }
        , measurement_id{ value_or(context.header.measurementInformation->measurementID, ""s) }, measurement_storage(context.storage.measurement) {

        if (!perform_noise_adjust)
            return;

        GDEBUG("Folder to store noise dependencies is %s\n", noise_dependency_folder.c_str());
        GDEBUG("NoiseAdjustGadget::perform_noise_adjust_ is %d\n", perform_noise_adjust);
        GDEBUG("NoiseAdjustGadget::pass_nonconformant_data_ is %d\n", pass_nonconformant_data);
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
            size_t CHA = noise_covariance->matrix_.get_size(0);
            if (noise_covariance->labels_.size() == CHA) {
                std::vector<std::string> current_coil_labels;
                if (current_ismrmrd_header.acquisitionSystemInformation) {
                    for (auto& l : current_ismrmrd_header.acquisitionSystemInformation->coilLabel) {
                        current_coil_labels.push_back(l.coilName);
                    }
                }

                std::vector<size_t> coil_order_of_data_in_noise;
                bool labels_match = compare_coil_label(noise_covariance->labels_,
                    current_coil_labels, coil_order_of_data_in_noise);

                if (!labels_match) {
                    // if number of channels in noise is different than data
                    // or
                    // if any channels in noise do not exist in data
                    if (CHA != current_coil_labels.size()) {
                        GDEBUG("Noise and measurement have different number of coils\n");
                    } else {
                        if (coil_order_of_data_in_noise.size() == CHA) {
                            GWARN_STREAM("Noise and meansurement have different coils, but will be reordered ... ");
                            noise_covariance->matrix_ = reorder_noise_channels(
                                noise_covariance->matrix_, coil_order_of_data_in_noise);

                        } else {
                            GWARN_STREAM("Noise and meansurement have different coils and cannot be reordered ... ");
                        }
                    }
                }
                loaded_noise_labels = noise_covariance->labels_;
                loaded_noise_sample_count = noise_covariance->sample_count_;
                GDEBUG("\033[32mSuccessfully loaded stored noise data with %zu channels and %zu samples\033[0m\n", CHA, noise_covariance->sample_count_);
                // print_covariance_matrix(noise_covariance->matrix_, loaded_noise_labels, loaded_noise_sample_count);

                return LoadedNoise{noise_covariance->matrix_, noise_covariance->noise_dwell_time_us_};

            } else if (current_ismrmrd_header.acquisitionSystemInformation) {
                GERROR("\033[31mNoise covariance matrix is malformed. Number of labels does not match number of channels.\033[0m\n");
            }
        }

        // No noise data found, gather it
        GDEBUG("\033[31mNo noise covariance data found, will gather from noise scans.\033[0m\n");
        return NoiseGatherer{};
    }

    static bool is_noise(const Core::Acquisition& acq) {
        return std::get<ISMRMRD::AcquisitionHeader>(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
    }

    template <class NOISEHANDLER>
    void NoiseAdjustGadget::add_noise(NOISEHANDLER& nh, const Gadgetron::Core::Acquisition&) const {
    }

    template <> void NoiseAdjustGadget::add_noise(NoiseGatherer& ng, const Gadgetron::Core::Acquisition& acq) const {
        auto& data    = std::get<hoNDArray<std::complex<float>>>(acq);
        auto& head    = std::get<ISMRMRD::AcquisitionHeader>(acq);
        if (ng.tmp_covariance.empty()) {
            auto channels = head.active_channels;
            ng.tmp_covariance = hoNDArray<std::complex<float>>(channels, channels);
            std::fill(ng.tmp_covariance.begin(), ng.tmp_covariance.end(), std::complex<float>(0));
        }

        if (ng.noise_dwell_time_us == 0)
            ng.noise_dwell_time_us = head.sample_time_us;

        auto dataM = as_arma_matrix(data);
        auto covariance = as_arma_matrix(ng.tmp_covariance);
        covariance += dataM.t()*dataM;


        ng.total_number_of_samples += head.number_of_samples;
    }

    template <> void NoiseAdjustGadget::add_noise(NoiseHandler& nh, const Gadgetron::Core::Acquisition& acq) const {
        Core::visit([&](auto& var) { this->add_noise(var, acq); }, nh);
    }
    template <class NOISEHANDLER> void NoiseAdjustGadget::save_noisedata(NOISEHANDLER& nh) {}

    template <> void NoiseAdjustGadget::save_noisedata(NoiseGatherer& ng) {
        if (ng.tmp_covariance.empty())
            return;

        normalize_covariance(ng);

        std::vector<std::string> coil_labels;
        for (auto& label : current_ismrmrd_header.acquisitionSystemInformation->coilLabel) {
            coil_labels.push_back(label.coilName);
        }

        auto noise_covariance = NoiseCovariance( 
            ng.tmp_covariance.get_size(0),
            coil_labels,
            ng.tmp_covariance,
            ng.total_number_of_samples,
            ng.noise_dwell_time_us,
            receiver_noise_bandwidth);

        if (!noise_covariance_out.empty()) {
            std::ofstream os(noise_covariance_out, std::ios::out | std::ios::binary);
            if (os.is_open()) {
                GDEBUG("\033[32mWriting noise covariance from %zu samples to %s\033[0m\n", ng.total_number_of_samples, noise_covariance_out.c_str());
                noise_covariance.SerializeToSfndam(os);
                os.flush();
                os.close();
            } else {
                GERROR("Unable to open file %s for writing noise covariance\n", noise_covariance_out.c_str());
            }
        } else {
            GDEBUG("\033[32mSaving noise covariance from %zu samples to storage server\033[0m\n", ng.total_number_of_samples);
            this->measurement_storage->store("noise_covariance", noise_covariance);
        }
    }

    template <> void NoiseAdjustGadget::save_noisedata(NoiseHandler& nh) {
        Core::visit([&](auto& var) { this->save_noisedata(var); }, nh);
    }


    template <class NH>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(NH nh, Core::Acquisition& acq) {
        return std::move(nh);
    };

    template <>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(
        Prewhitener pw, Core::Acquisition& acq) {

        auto& data = std::get<hoNDArray<std::complex<float>>>(acq);
        if (data.get_size(1) == pw.prewhitening_matrix.get_size(0)) {
            auto dataM = as_arma_matrix(data);
            auto pwm = as_arma_matrix(pw.prewhitening_matrix);
            dataM *= pwm;
        } else if (!this->pass_nonconformant_data) {
            throw std::runtime_error("Input data has different number of channels from noise data");
        }
        return std::move(pw);
    }

    template <>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(
        NoiseGatherer ng, Core::Acquisition& acq) {
        auto& head = std::get<ISMRMRD::AcquisitionHeader>(acq);
        if (ng.total_number_of_samples == 0)
            return std::move(ng);

        this->save_noisedata(ng);

        auto masked_covariance = mask_channels(ng.tmp_covariance, scale_only_channels);

        auto prewhitening_matrix = computeNoisePrewhitener(masked_covariance);
        prewhitening_matrix
            *= calculate_scale_factor(head.sample_time_us, ng.noise_dwell_time_us, receiver_noise_bandwidth);
        return handle_acquisition(Prewhitener{ prewhitening_matrix }, acq);
    }

    template <>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(
        LoadedNoise ln, Core::Acquisition& acq)  {
        auto& head               = std::get<ISMRMRD::AcquisitionHeader>(acq);
        auto masked_covariance   = mask_channels(std::move(ln.covariance), scale_only_channels);
        auto prewhitening_matrix = computeNoisePrewhitener(masked_covariance);
        prewhitening_matrix
            *= calculate_scale_factor(head.sample_time_us, ln.noise_dwell_time_us, receiver_noise_bandwidth);
        return handle_acquisition(Prewhitener{ prewhitening_matrix }, acq);
    }

    template <>
    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::handle_acquisition(
        NoiseHandler nh, Core::Acquisition& acq) {
        return Core::visit([&](auto var) { return this->handle_acquisition<decltype(var)>(std::move(var), acq); }, std::move(nh));
    }

    void NoiseAdjustGadget::process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) {

        scale_only_channels = current_ismrmrd_header.acquisitionSystemInformation
                                  ? find_scale_only_channels(scale_only_channels_by_name,
                                      current_ismrmrd_header.acquisitionSystemInformation->coilLabel)
                                  : std::vector<size_t>{};


        bool first_non_noise = true;
        for (auto acq : input) {
            if (is_noise(acq)) {
                add_noise(noisehandler, acq);
                continue;
            }
            if (first_non_noise) {
                first_non_noise = false;
                Core::visit([this](const auto& h) {
                    if constexpr (std::is_same_v<std::decay_t<decltype(h)>, NoiseGatherer>) {
                        GDEBUG("Processing first non-noise scan.  Using noise scans gathered during this scan\n");
                        std::vector<std::string> coil_labels;
                        if (current_ismrmrd_header.acquisitionSystemInformation)
                            for (auto& l : current_ismrmrd_header.acquisitionSystemInformation->coilLabel)
                                coil_labels.push_back(l.coilName);
                        print_covariance_matrix(h.tmp_covariance, coil_labels, h.total_number_of_samples);
                    } else if constexpr (std::is_same_v<std::decay_t<decltype(h)>, LoadedNoise>) {
                        GDEBUG("Processing first non-noise scan.  Using loaded covariance matrix from noise scans gathered in a dependent measurement\n");
                        print_covariance_matrix(h.covariance, loaded_noise_labels, loaded_noise_sample_count);
                    }
                }, noisehandler);
            }
            noisehandler = handle_acquisition(std::move(noisehandler), acq);
            output.push(std::move(acq));
        }

        this->save_noisedata(noisehandler);
    }

    Core::optional<NoiseCovariance> NoiseAdjustGadget::load_noisedata() const {
        if (!noise_covariance_in.empty()) {
            GDEBUG("\033[33mAttempting to load noise covariance file %s\033[0m\n", noise_covariance_in.c_str());
            std::ifstream file(noise_covariance_in, std::ios::binary);
            if (!file) {
                GERROR("Could not open noise covariance file %s\n", noise_covariance_in.c_str());
                throw std::runtime_error("Could not open noise covariance file");
            }
            return NoiseCovariance::DeserializeFromSfnadm(file);
        } else {
            GDEBUG("Measurement ID is %s\n", measurement_id.c_str());
            if (!current_ismrmrd_header.measurementInformation) {
                GWARN("ISMRMRD Header is missing measurmentinformation. Skipping noise adjust");
                return Core::none;
            }
            const auto& measurementDependency = current_ismrmrd_header.measurementInformation->measurementDependency;
            auto val = std::find_if(measurementDependency.begin(), measurementDependency.end(), [](const auto& dependency) {
                return boost::algorithm::to_lower_copy(dependency.dependencyType) == "noise";
            });

            if (val == measurementDependency.end())
                return Core::none;

            auto noise_dependency = *val;
            GDEBUG("\033[33mAttempting to retrieve noise_covariance from storage server, measurementID %s\033[0m\n", noise_dependency.measurementID.c_str());
            return measurement_storage->get_latest<NoiseCovariance>(noise_dependency.measurementID, "noise_covariance");
        }
    }

    GADGETRON_GADGET_EXPORT(NoiseAdjustGadget)

} // namespace Gadgetron
