
#include "CmrParametricT1SRMappingGadget.h"
#include <iomanip>
#include <sstream>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_utility.h"
#include "cmr_t1_mapping.h"

namespace Gadgetron {

    CmrParametricT1SRMappingGadget::CmrParametricT1SRMappingGadget() : BaseClass()
    {
    }

    CmrParametricT1SRMappingGadget::~CmrParametricT1SRMappingGadget()
    {
    }

    int CmrParametricT1SRMappingGadget::process_config(const mrd::Header& header)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(header) == GADGET_OK, GADGET_FAIL);

        auto& h = header;

        if (!h.acquisition_system_information)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        if (this->imaging_prep_time_from_protocol.value())
        {
            this->prep_times_ = h.sequence_parameters->t_i; // TI is in the unit of seconds
            GDEBUG_STREAM("Read prep times from from protocol : " << this->prep_times_.size() << " [ ");
            for (size_t n = 0; n < this->prep_times_.size(); n++)
            {
                this->prep_times_[n] *= 1000; // convert to ms
                GDEBUG_STREAM(this->prep_times_[n]);
            }
            GDEBUG_STREAM(" ] ");
        }

        // -------------------------------------------------

        return GADGET_OK;
    }

    int CmrParametricT1SRMappingGadget::perform_mapping(mrd::ImageArray& data, mrd::ImageArray& map, mrd::ImageArray& para, mrd::ImageArray& map_sd, mrd::ImageArray& para_sd)
    {
        try
        {
            if (perform_timing.value()) { gt_timer_.start("CmrParametricT1SRMappingGadget::perform_mapping"); }

            GDEBUG_CONDITION_STREAM(verbose.value(), "CmrParametricT1SRMappingGadget::perform_mapping(...) starts ... ");

            size_t RO = data.data.get_size(0);
            size_t E1 = data.data.get_size(1);
            size_t E2 = data.data.get_size(2);
            size_t CHA = data.data.get_size(3);
            size_t N = data.data.get_size(4);
            size_t S = data.data.get_size(5);
            size_t SLC = data.data.get_size(6);

            size_t ro, e1, s, slc, p;

            GADGET_CHECK_RETURN(E2 == 1, GADGET_FAIL);
            GADGET_CHECK_RETURN(CHA == 1, GADGET_FAIL);
            GADGET_CHECK_RETURN(this->prep_times_.size() >= N, GADGET_FAIL);

            hoNDArray<float> mag;
            Gadgetron::abs(data.data, mag);

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array(mag, debug_folder_full_path_ + "CmrParametricT1SRMapping_data_mag");
            }

            bool need_sd_map = send_sd_map.value();

            Gadgetron::GadgetronTimer gt_timer(false);

            // -------------------------------------------------------------
            // set mapping parameters

            Gadgetron::CmrT1SRMapping<float> t1_sr;

            t1_sr.fill_holes_in_maps_ = perform_hole_filling.value();
            t1_sr.max_size_of_holes_ = max_size_hole.value();
            t1_sr.compute_SD_maps_ = need_sd_map;

            t1_sr.ti_.resize(N, 0);
            memcpy(&(t1_sr.ti_)[0], &this->prep_times_[0], sizeof(float)*N);

            // set the anchor image TS
            size_t anchor_ind = this->anchor_image_index.value();
            if (anchor_ind < N)
            {
                t1_sr.ti_[anchor_ind] = this->anchor_TS.value();
            }

            t1_sr.data_.create(RO, E1, N, S, SLC, mag.begin());

            t1_sr.max_iter_ = max_iter.value();
            t1_sr.thres_fun_ = thres_func.value();
            t1_sr.max_map_value_ = max_T1.value();

            t1_sr.verbose_ = verbose.value();
            t1_sr.debug_folder_ = debug_folder_full_path_;
            t1_sr.perform_timing_ = perform_timing.value();

            // -------------------------------------------------------------
            // compute mask if needed
            if (mapping_with_masking.value())
            {
                t1_sr.mask_for_mapping_.create(RO, E1, SLC);

                // get the image with longest TS time
                hoNDArray<float> mag_longest_TS;
                mag_longest_TS.create(RO, E1, SLC);

                for (slc = 0; slc < SLC; slc++)
                {
                    size_t ind = N - 1;
                    if (anchor_ind < N)
                    {
                        ind = anchor_ind;
                    }
                    else
                    {
                        float max_ts = this->prep_times_[0];
                        for (size_t n = 1; n < this->prep_times_.size(); n++)
                        {
                            if(this->prep_times_[n]>max_ts)
                            {
                                max_ts = this->prep_times_[n];
                                ind = n;
                            }
                        }
                    }

                    memcpy(&mag_longest_TS(0, 0, slc), &mag(0, 0, ind, 0, slc), sizeof(float)*RO*E1);
                }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(mag_longest_TS, debug_folder_full_path_ + "CmrParametricT1SRMapping_mag_longest_TS");
                }

                double scale_factor = 1.0;
                if (data.meta[0].count(GADGETRON_IMAGE_SCALE_RATIO) > 0)
                {
                    scale_factor = std::get<double>(data.meta[0][GADGETRON_IMAGE_SCALE_RATIO].front());
                }

                GDEBUG_STREAM("CmrParametricT1SRMappingGadget, find incoming image has scale factor of " << scale_factor);

                if (perform_timing.value()) { gt_timer.start("CmrParametricT1SRMappingGadget::compute_mask_for_mapping"); }
                this->compute_mask_for_mapping(mag, t1_sr.mask_for_mapping_, (float)scale_factor);
                if (perform_timing.value()) { gt_timer.stop(); }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(t1_sr.mask_for_mapping_, debug_folder_full_path_ + "CmrParametricT1SRMapping_mask_for_mapping");
                }
            }

            // -------------------------------------------------------------
            // perform mapping

            if (perform_timing.value()) { gt_timer.start("CmrParametricT1SRMappingGadget, t1_sr.perform_parametric_mapping"); }
            t1_sr.perform_parametric_mapping();
            if (perform_timing.value()) { gt_timer.stop(); }

            size_t num_para = t1_sr.get_num_of_paras();

            // -------------------------------------------------------------
            // get the results

            map.data.create(RO, E1, E2, CHA, 1, S, SLC);
            Gadgetron::clear(map.data);
            map.headers.create(1, S, SLC);
            map.meta.create(1, S, SLC);

            para.data.create(RO, E1, E2, CHA, num_para, S, SLC);
            Gadgetron::clear(para.data);
            para.headers.create(num_para, S, SLC);
            para.meta.create(num_para, S, SLC);

            if (need_sd_map)
            {
                map_sd.data.create(RO, E1, E2, CHA, 1, S, SLC);
                Gadgetron::clear(map_sd.data);
                map_sd.headers.create(1, S, SLC);
                map_sd.meta.create(1, S, SLC);

                para_sd.data.create(RO, E1, E2, CHA, num_para, S, SLC);
                Gadgetron::clear(para_sd.data);
                para_sd.headers.create(num_para, S, SLC);
                para_sd.meta.create(num_para, S, SLC);
            }

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (e1 = 0; e1 < E1; e1++)
                    {
                        for (ro = 0; ro < RO; ro++)
                        {
                            map.data(ro, e1, 0, 0, 0, s, slc) = t1_sr.map_(ro, e1, s, slc);

                            if (need_sd_map)
                            {
                                map_sd.data(ro, e1, 0, 0, 0, s, slc) = t1_sr.sd_map_(ro, e1, s, slc);
                            }

                            for (p = 0; p < num_para; p++)
                            {
                                para.data(ro, e1, 0, 0, p, s, slc) = t1_sr.para_(ro, e1, p, s, slc);

                                if (need_sd_map)
                                {
                                    para_sd.data(ro, e1, 0, 0, p, s, slc) = t1_sr.sd_para_(ro, e1, p, s, slc);
                                }
                            }
                        }
                    }

                    size_t slc_ind = data.headers(0, s, slc).slice.value_or(0);

                    map.headers(0, s, slc) = data.headers(0, s, slc);
                    map.headers(0, s, slc).image_index = 1 + slc_ind;
                    map.headers(0, s, slc).image_series_index = 11;
                    map.meta(0, s, slc) = data.meta(0, s, slc);

                    auto& map_meta = map.meta(0, s, slc);
                    map_meta[GADGETRON_DATA_ROLE] = { GADGETRON_IMAGE_T1MAP };
                    map_meta[GADGETRON_SEQUENCEDESCRIPTION].push_back(GADGETRON_IMAGE_T1MAP);
                    map_meta[GADGETRON_IMAGEPROCESSINGHISTORY].push_back(GADGETRON_IMAGE_T1MAP);

                    map_sd.headers(0, s, slc) = data.headers(0, s, slc);
                    map_sd.headers(0, s, slc).image_index = 1 + slc_ind;
                    map_sd.headers(0, s, slc).image_series_index = 12;
                    map_sd.meta(0, s, slc) = data.meta(0, s, slc);

                    auto& map_sd_meta = map_sd.meta(0, s, slc);
                    map_sd_meta[GADGETRON_DATA_ROLE] = { GADGETRON_IMAGE_T1SDMAP };
                    map_sd_meta[GADGETRON_SEQUENCEDESCRIPTION].push_back(GADGETRON_IMAGE_T1SDMAP);
                    map_sd_meta[GADGETRON_IMAGEPROCESSINGHISTORY].push_back(GADGETRON_IMAGE_T1SDMAP);

                    if (need_sd_map)
                    {
                        for (p = 0; p < num_para; p++)
                        {
                            para.headers(p, s, slc) = data.headers(0, s, slc);
                            para.headers(p, s, slc).image_index = 1 + p + slc_ind*num_para;
                            para.meta(p, s, slc) = data.meta(0, s, slc);

                            para_sd.headers(p, s, slc) = data.headers(0, s, slc);
                            para_sd.headers(p, s, slc).image_index = 1 + p + slc_ind*num_para;
                            para_sd.meta(p, s , slc) = data.meta(0, s, slc);
                        }
                    }
                }
            }

            // -------------------------------------------------------------

            if (perform_timing.value()) { gt_timer_.stop(); }
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in CmrParametricT1SRMappingGadget::perform_mapping(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int CmrParametricT1SRMappingGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "CmrParametricT1SRMappingGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(CmrParametricT1SRMappingGadget)

}
