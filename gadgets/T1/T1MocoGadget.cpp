//
// Created by dchansen on 3/9/20.
//
#include "demons_registration.h"
#include <unordered_set>

#include "PureGadget.h"
#include "cmr_parametric_mapping.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"
#include "mri_core_data.h"
#include "mri_core_def.h"
#include "t1fit.h"
#include <range/v3/algorithm.hpp>
#include <range/v3/range_concepts.hpp>
namespace Gadgetron {

class T1MocoGadget : public Core::ChannelGadget<IsmrmrdImageArray> {

  public:
    T1MocoGadget(const Core::Context& context, const Core::GadgetProperties& properties)
        : Core::ChannelGadget<IsmrmrdImageArray>(context, properties),
          field_strength{*context.header.acquisitionSystemInformation->systemFieldStrength_T} {}

    NODE_PROPERTY(correction_factor, float, "Empirical correction factor for T1", 1.0365f);
    NODE_PROPERTY(regularization_sigma, float, "Gaussian regularizatin for registration", 2.0f);
    NODE_PROPERTY(demons_iterations, unsigned int, "Number of iterations for the demons registration", 40);
    NODE_PROPERTY(step_size, float, "Maximum step size for demons registration (between 0.1 and 2.0)", 2.0f);
    NODE_PROPERTY(iterations, unsigned int, "Number of iterations of demons registration and T1 fit", 5);
    NODE_PROPERTY(scales, unsigned int, "Number of image scales to use", 1);

  private:
    void process(Core::InputChannel<IsmrmrdImageArray>& input, Core::OutputChannel& out) final override {

        GDEBUG("Sigma %f demons_iterations %i Step size %f iterations %i\n", regularization_sigma, demons_iterations,
               step_size, iterations);

        for (auto images : input) {

            auto TI_values = extract_MOLLI_TI(*images.acq_headers_);
            auto data_dims = images.data_.dimensions();
            images.data_.reshape(data_dims[0], data_dims[1], -1);

            sort_images_and_values(images, TI_values);

            auto moco_images = multi_stage_T1_registration(images.data_, TI_values);

            auto phase_corrected = T1::phase_correct(moco_images, TI_values);

            const auto [A, B, T1star] = T1::fit_T1_3param(phase_corrected, TI_values);

            auto T1 = t1_from_t1star(A, B, T1star);
            clean_image(T1);
            perform_hole_filling(T1);

            auto error_map = T1::calculate_error_map({A, B, T1}, phase_corrected, TI_values);
            clean_image(error_map,1000);

            T1 *= correction_factor;

            auto header = images.headers_[0];
            header.data_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
            header.image_series_index = 5;
            auto meta = create_T1_meta(images.meta_.front(),T1);
            auto sd_meta = create_T1SD_meta(images.meta_.front());
            // send original images
            images.data_.reshape(data_dims);
            set_RAW_headers_and_meta(images, TI_values);
            out.push(images);

            // send MOCO images
            images.data_ = hoNDArray<std::complex<float>>(phase_corrected);
            images.data_.reshape(data_dims);

            auto mag_images = images;
            mag_images.data_ = hoNDArray<std::complex<float>>(abs(moco_images));
            mag_images.data_.reshape(data_dims);
            set_MOCO_MAG_headers_and_meta(mag_images, TI_values);
            out.push(std::move(mag_images));

            set_PSIR_headers_and_meta(images, TI_values);
            out.push(std::move(images));

            // send out T1 map
            auto sd_header = header;
            sd_header.image_series_index = 4;
            out.push(Core::Image<float>{sd_header, std::move(error_map), sd_meta});
            out.push(Core::Image<float>{header, std::move(T1), meta});
        }
    }

    ISMRMRD::MetaContainer create_T1SD_meta(ISMRMRD::MetaContainer meta) const {

        double scaling_factor = 1;
        double window_center = 200;
        double window_width = 400;
        std::string lut =
            std::abs(field_strength - float(1.5)) < 1e-1 ? "GadgetronT1_IR_1_5T.pal" : "GadgetronT1_IR_3T.pal";

        std::ostringstream ostr;
        ostr << "x" << scaling_factor;
        std::string scalingStr = ostr.str();

        std::ostringstream ostr_unit;
        ostr_unit << std::setprecision(3) << 1.0f / scaling_factor << "ms";
        std::string unitStr = ostr_unit.str();

        meta.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1SDMAP);
        meta.append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1SDMAP);
        meta.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_MOCO);
        meta.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1SDMAP);

        meta.set(GADGETRON_IMAGE_SCALE_RATIO, scaling_factor);
        meta.set(GADGETRON_IMAGE_WINDOWCENTER, (long)(window_center * scaling_factor));
        meta.set(GADGETRON_IMAGE_WINDOWWIDTH, (long)(window_width * scaling_factor));
        meta.set(GADGETRON_IMAGE_COLORMAP, lut.c_str());

        meta.set(GADGETRON_IMAGECOMMENT, meta.as_str(GADGETRON_DATA_ROLE));
        meta.append(GADGETRON_IMAGECOMMENT, scalingStr.c_str());
        meta.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());

        return std::move(meta);
    }

    static std::pair<double,double> find_window_center_width_T1(const hoNDArray<float>& image, const ISMRMRD::IsmrmrdHeader& header ){

        std::pair<double,double> pre_window = {1300,1300};
        std::pair<double, double> post_window = {400,300};

        try {
            auto protocolname = header.measurementInformation.get().protocolName.get();

            std::regex is_post(R"((?:^|[^A-Z])post(?:$|[^A-Z]))", std::regex_constants::icase | std::regex_constants::ECMAScript);
            if (std::regex_search(protocolname,is_post)) return post_window;

            std::regex is_pre(R"((?:^|[^A-Z])pre(?:$|[^A-Z]))", std::regex_constants::icase | std::regex_constants::ECMAScript);
            if (std::regex_search(protocolname,is_pre)) return pre_window;

        } catch (...) {}


        auto top = percentile(image, 0.9f);

        auto window_distance = [top](const auto& window){ return (top - ( window.first+window.second/2));};

        if (window_distance(pre_window) < window_distance(post_window))
            return pre_window;
        else
            return post_window;
    }

    ISMRMRD::MetaContainer create_T1_meta(ISMRMRD::MetaContainer meta, const hoNDArray<float>& t1map) const {

        double scaling_factor = 1;

        auto [window_center, window_width] = find_window_center_width_T1(t1map,header);
        GDEBUG("Setting T1 window level to %f %f \n", window_center, window_width);

        std::string lut =
            std::abs(field_strength - float(1.5)) < 1e-1 ? "GadgetronT1_IR_1_5T.pal" : "GadgetronT1_IR_3T.pal";

        std::ostringstream ostr;
        ostr << "x" << scaling_factor;
        std::string scalingStr = ostr.str();

        std::ostringstream ostr_unit;
        ostr_unit << std::setprecision(3) << 1.0f / scaling_factor << "ms";
        std::string unitStr = ostr_unit.str();

        meta.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1MAP);
        meta.append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1MAP);
        meta.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_MOCO);
        meta.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1MAP);

        meta.set(GADGETRON_IMAGE_SCALE_RATIO, scaling_factor);
        meta.set(GADGETRON_IMAGE_WINDOWCENTER, (long)(window_center * scaling_factor));
        meta.set(GADGETRON_IMAGE_WINDOWWIDTH, (long)(window_width * scaling_factor));
        meta.set(GADGETRON_IMAGE_COLORMAP, lut.c_str());

        meta.set(GADGETRON_IMAGECOMMENT, meta.as_str(GADGETRON_DATA_ROLE));
        meta.append(GADGETRON_IMAGECOMMENT, scalingStr.c_str());
        meta.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());

        return std::move(meta);
    }

    static std::vector<float> extract_MOLLI_TI(const hoNDArray<ISMRMRD::AcquisitionHeader>& acq_headers) {

        std::map<int, std::vector<ISMRMRD::AcquisitionHeader>> look_locker_sets;

        for (auto subarray : spans(acq_headers, 1)) {
            auto& header = *std::find_if(subarray.begin(), subarray.end(), [](ISMRMRD::AcquisitionHeader& acq_header) {
                return acq_header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
            });
            look_locker_sets[header.user_int[5]].push_back(header);
        }

        std::vector<float> TI_values;

        for (auto& [set, headers] : look_locker_sets) {
            auto minimum_acquisition_time_stamp = std::accumulate(
                headers.begin(), headers.end(), std::numeric_limits<uint32_t>::max(), [](auto val1, auto val2) {
                    return val1 < val2.acquisition_time_stamp ? val1 : val2.acquisition_time_stamp;
                });

            for (auto& header : headers) {
                float TI_value =
                    (header.acquisition_time_stamp - minimum_acquisition_time_stamp) * 2.5f + header.user_int[4];
                TI_values.push_back(TI_value);
                GDEBUG("set %d look-locker %d ti: %f  acq_time_stamp: %d  \n", header.idx.set, set, TI_value,
                       header.acquisition_time_stamp);
            }
        }

        return TI_values;
    }

    static void clean_image(hoNDArray<float>& data, float upper_limit = 5000.0f) {
        std::transform(data.begin(), data.end(), data.begin(), [upper_limit](auto val) {
            if (val <= 0)
                return 0.0f;
            if (val >= upper_limit)
                return 0.0f;
            if (std::isnan(val))
                return 0.0f;
            return val;
        });
    }

    void sort_images_and_values(IsmrmrdImageArray& images, std::vector<float>& TI_values) {
        auto sorted_index = argsort(TI_values);

        auto dims = images.data_.dimensions();
        images.data_.reshape(dims[0], dims[1], -1);
        auto data_copy = images.data_;
        std::fill(images.data_.begin(), images.data_.end(), std::complex<float>(1));

        for (size_t i = 0; i < data_copy.get_size(2); i++) {
            using namespace Gadgetron::Indexing;
            images.data_(slice, slice, i) = data_copy(slice, slice, sorted_index[i]);
        }
        images.data_.reshape(dims);

        auto TI_sorted = TI_values;

        for (size_t i = 0; i < TI_values.size(); i++)
            TI_sorted[i] = TI_values[sorted_index[i]];

        TI_values = std::move(TI_sorted);
    }
    static void set_RAW_headers_and_meta(IsmrmrdImageArray& images, const std::vector<float>& TI_values) {
        for (auto& header : images.headers_) {
            header.image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
            header.image_series_index = 1;
        }

        for (size_t ind = 0; ind < images.meta_.size(); ind++) {
            auto& meta = images.meta_[ind];
            meta.set(GADGETRON_IMAGE_INVERSIONTIME, (double)(TI_values[ind]));
        }
    }

    static void set_MOCO_MAG_headers_and_meta(IsmrmrdImageArray& images, const std::vector<float>& TI_values) {
        for (auto& header : images.headers_) {
            header.image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
            header.image_series_index = 2;
        }

        for (size_t ind = 0; ind < images.meta_.size(); ind++) {
            auto& meta = images.meta_[ind];
            meta.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_MOCO);
            meta.set(GADGETRON_IMAGE_INVERSIONTIME, (double)(TI_values[ind]));
            meta.append(GADGETRON_SEQUENCEDESCRIPTION,GADGETRON_IMAGE_MOCO);
        }
    }

    static void set_PSIR_headers_and_meta(IsmrmrdImageArray& images, const std::vector<float>& TI_values) {
        for (auto& header : images.headers_) {
            header.image_type = ISMRMRD::ISMRMRD_IMTYPE_REAL;
            header.image_series_index = 3;
        }

        images.data_ += 4096.0;
        for (size_t ind = 0; ind < images.meta_.size(); ind++) {
            auto& meta = images.meta_[ind];
            meta.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_PSIR);
            meta.append(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_MOCO);
            meta.set(GADGETRON_IMAGE_SCALE_OFFSET, (double)(4096.0));
            meta.set(GADGETRON_IMAGE_INVERSIONTIME, (double)(TI_values[ind]));
            meta.append(GADGETRON_SEQUENCEDESCRIPTION,GADGETRON_IMAGE_MOCO);
            meta.append(GADGETRON_SEQUENCEDESCRIPTION,GADGETRON_IMAGE_PSIR);
        }
    }

    template <class T> std::vector<size_t> argsort(const std::vector<T>& data) {

        std::vector<std::pair<size_t, T>> index_and_values(data.size());

        for (size_t i = 0; i < data.size(); i++)
            index_and_values[i] = {i, data[i]};

        std::stable_sort(index_and_values.begin(), index_and_values.end(),
                         [](auto val1, auto val2) { return val1.second < val2.second; });

        std::vector<size_t> index(data.size());

        for (size_t i = 0; i < data.size(); i++)
            index[i] = index_and_values[i].first;
        return index;
    }

    hoNDArray<std::complex<float>> multi_stage_T1_registration(const hoNDArray<std::complex<float>>& data,
                                                               const std::vector<float>& TIs) const {

        auto abs_data = abs(data);
        auto first_vfields = register_compatible_frames(abs_data, TIs);
        auto second_vfields = T1::t1_registration(data, TIs, std::move(first_vfields), iterations,
                                                  {demons_iterations, regularization_sigma, step_size});

        auto deformed_images = T1::deform_groups(abs_data, second_vfields);
        auto final_vfield = T1::register_groups_CMR(abs_data, deformed_images, 24.0f);
        return T1::deform_groups(data, final_vfield);
    }

    hoNDArray<vector_td<float, 2>> register_compatible_frames(const hoNDArray<float>& abs_data,
                                                              const std::vector<float>& TIs) const {
        using namespace Gadgetron::Indexing;
        using namespace ranges;
        auto arg_max_TI = max_element(TIs) - TIs.begin();

        const hoNDArray<float> reference_frame = abs_data(slice, slice, arg_max_TI);

        const auto valid_transforms =
            view::iota(size_t(0), abs_data.get_size(2)) |
            views::filter([&arg_max_TI](auto index) { return index != arg_max_TI; }) | views::filter([&](auto index) {
                return jensen_shannon_divergence(abs_data(slice, slice, index), reference_frame) < 0.2;
            }) |
            to<std::vector>();

        std::vector<hoNDArray<vector_td<float, 2>>> vfields(abs_data.get_size(2));

#pragma omp parallel for default(shared)
        for (int i = 0; i < int(valid_transforms.size()); i++) {
            vfields[valid_transforms[i]] =
                T1::register_groups_CMR(abs_data(slice, slice, valid_transforms[i]), reference_frame);
        }
        auto missing_indices = view::iota(size_t(0), abs_data.get_size(2)) | view::filter([&](auto index) {
                                   return !binary_search(valid_transforms, index) && (index != arg_max_TI);
                               });

        for (auto index : missing_indices) {
            auto closest_index = lower_bound(valid_transforms, index);
            if (closest_index != valid_transforms.end()) {
                vfields[index] = vfields[*closest_index];
            } else {
                vfields[index] = hoNDArray<vector_td<float, 2>>(abs_data.get_size(0), abs_data.get_size(1), 1);
                vfields[index].fill(vector_td<float, 2>(0, 0));
            }
        }

        vfields[arg_max_TI] = hoNDArray<vector_td<float, 2>>(abs_data.get_size(0), abs_data.get_size(1), 1);
        vfields[arg_max_TI].fill(vector_td<float, 2>(0, 0));
        return concat_along_dimension(vfields, 2);
    }

    hoNDArray<float> t1_from_t1star(const hoNDArray<float>& A, const hoNDArray<float>& B,
                                    const hoNDArray<float>& T1star) {

        auto T1 = hoNDArray<float>(T1star.dimensions());

        for (size_t i = 0; i < T1.size(); i++)
            T1[i] = T1star[i] * (B[i] / A[i] - 1.0f);

        return T1;
    }

    float field_strength;
};

GADGETRON_GADGET_EXPORT(T1MocoGadget)
} // namespace Gadgetron