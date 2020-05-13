//
// Created by dchansen on 3/9/20.
//
#include "demons_registration.h"
#include <unordered_set>

#include "PureGadget.h"
#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_fileio.h"
#include "mri_core_data.h"
#include "mri_core_def.h"
#include "cmr_parametric_mapping.h"
#include "t1fit.h"

namespace Gadgetron {

    class T1Gadget : public Core::PureGadget<Core::Image<float>, IsmrmrdImageArray> {

    public:
        T1Gadget(const Core::Context& context, const Core::GadgetProperties& properties)
            : Core::PureGadget<Core::Image<float>, IsmrmrdImageArray>(context, properties)
            , TIs{ *(context.header.sequenceParameters->TI) } { }

        NODE_PROPERTY(correction_factor, float, "Empirical correction factor for T1", 1.0365f);

    private:
        std::vector<float> extract_MOLLI_TI(const hoNDArray<ISMRMRD::AcquisitionHeader>& acq_headers) const {

            std::map<int, std::vector<ISMRMRD::AcquisitionHeader>> look_locker_sets;

            for (auto subarray : spans(acq_headers, 1)) {
                auto& header
                    = *std::find_if(subarray.begin(), subarray.end(), [](ISMRMRD::AcquisitionHeader& acq_header) {
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
                    float TI_value = (header.acquisition_time_stamp - minimum_acquisition_time_stamp) * 2.5f
                                    + header.user_int[4];
                    TI_values.push_back(TI_value);
                    GDEBUG("set %d look-locker %d ti: %f  acq_time_stamp: %d  \n",header.idx.set,set,TI_value,header.acquisition_time_stamp);
                }
            }

            return TI_values;
        }

        static void clean_image(hoNDArray<float>& data){
            std::transform(data.begin(),data.end(),data.begin(),[](auto val){
                if (val <= 0) return 0.0f;
                if (val >= 5000) return 0.0f;
                if (std::isnan(val)) return 0.0f;
                return val;
            });
        }

        Core::Image<float> process_function(IsmrmrdImageArray images) const final override {

            auto TI_values = extract_MOLLI_TI(*images.acq_headers_);

            auto data_dims = images.data_.dimensions();
            images.data_.reshape({data_dims[0],data_dims[1],-1});
            auto vector_field = T1::t1_registration(images.data_, TI_values);

            auto moco_images = T1::deform_groups(images.data_, vector_field);

            auto phase_corrected = T1::phase_correct(moco_images, TI_values);


            auto [A, B, T1star] = T1::fit_T1_3param(phase_corrected, TI_values);

            B /= A;
            B -= 1;

            auto T1 = T1star;
            T1 *= B;
            T1 *= correction_factor;

            clean_image(T1);

            perform_hole_filling(T1);
            auto header               = images.headers_[0];
            header.image_series_index = ISMRMRD::ISMRMRD_IMTYPE_REAL;

            auto meta = images.meta_.front();
            meta.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1MAP);
            meta.append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1MAP);
            meta.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1MAP);
            return Core::Image<float>{ header, T1, Core::none };
        }

        const std::vector<float> TIs;
    };

    GADGETRON_GADGET_EXPORT(T1Gadget)
}