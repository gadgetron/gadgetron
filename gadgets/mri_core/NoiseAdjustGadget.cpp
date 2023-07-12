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

        // find the measurementID of this scan

        noisehandler = load_or_gather();

        if (context.parameters.find("noisecovariance") != context.parameters.end()) {
            noise_covariance_file_name = context.parameters.at("noisecovariance");
            GDEBUG_STREAM("Noise covariance matrix is provided as a parameter: " << noise_covariance_file_name);
        }
    }

    NoiseAdjustGadget::NoiseHandler NoiseAdjustGadget::load_or_gather() const {
        GDEBUG("Measurement ID is %s\n", measurement_id.c_str());
        if (!current_ismrmrd_header.measurementInformation) {
            GWARN("ISMRMRD Header is missing measurmentinformation. Skipping noise adjust");
            return NoiseGatherer{};
        }
        const auto& measurementDependency = current_ismrmrd_header.measurementInformation->measurementDependency;
        auto val = std::find_if(measurementDependency.begin(), measurementDependency.end(), [](const auto& dependency) {
            return boost::algorithm::to_lower_copy(dependency.dependencyType) == "noise";
        });

        // find the noise dependencies if any
        if (val == measurementDependency.end())
            return NoiseGatherer{};

        auto noise_dependency = *val;
        GDEBUG("Measurement ID of noise dependency is %s\n", noise_dependency.measurementID.c_str());

        auto noise_covariance = load_noisedata(noise_dependency.measurementID);

        // try to load the precomputed noise prewhitener
        if (!noise_covariance) {
            GDEBUG("Stored noise dependency is NOT found : %s\n", noise_dependency.measurementID.c_str());
            return NoiseGatherer{};
        } else {
            GDEBUG("Stored noise dependency is found : %s\n", noise_dependency.measurementID.c_str());
            GDEBUG("Stored noise dwell time in us is %f\n", noise_covariance->noise_dwell_time_us_);
            size_t CHA = noise_covariance->matrix_.get_size(0);
            GDEBUG("Stored noise channel number is %d\n", CHA);

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
                return LoadedNoise{noise_covariance->matrix_,noise_covariance->noise_dwell_time_us_};

            } else if (current_ismrmrd_header.acquisitionSystemInformation) {
                GERROR("Noise ismrmrd header does not have acquisition system information but current header "
                       "does\n");
            }

            //                    number_of_noise_samples_ = 1; // When we load the matrix, it is already
            //                    scaled.
        }
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

        this->measurement_storage->store("noise_covariance", noise_covariance);
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

    Core::optional<NoiseCovariance> NoiseAdjustGadget::load_noisedata(const std::string &noise_measurement_id) const {
       return measurement_storage->get_latest<NoiseCovariance>(noise_measurement_id, "noise_covariance");
    }

    GADGETRON_GADGET_EXPORT(NoiseAdjustGadget)

} // namespace Gadgetron
