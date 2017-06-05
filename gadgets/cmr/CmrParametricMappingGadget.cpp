
#include "CmrParametricMappingGadget.h"
#include <iomanip>
#include <sstream>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_utility.h"
#include "morphology.h"

namespace Gadgetron {

    CmrParametricMappingGadget::CmrParametricMappingGadget() : BaseClass()
    {
    }

    CmrParametricMappingGadget::~CmrParametricMappingGadget()
    {
    }

    int CmrParametricMappingGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!h.acquisitionSystemInformation)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        // -------------------------------------------------

        float field_strength_T_ = h.acquisitionSystemInformation.get().systemFieldStrength_T();
        GDEBUG_CONDITION_STREAM(verbose.value(), "field_strength_T_ is read from protocol : " << field_strength_T_);

        return GADGET_OK;
    }

    int CmrParametricMappingGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("CmrParametricMappingGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "CmrParametricMappingGadget::process(...) starts ... ");

        // -------------------------------------------------------------

        process_called_times_++;

        // -------------------------------------------------------------

        IsmrmrdImageArray* data = m1->getObjectPtr();

        // print out data info
        if (verbose.value())
        {
            GDEBUG_STREAM("----> CmrParametricMappingGadget::process(...) has been called " << process_called_times_ << " times ...");
            std::stringstream os;
            data->data_.print(os);
            GDEBUG_STREAM(os.str());
        }

        // -------------------------------------------------------------

        // some images do not need mapping
        if (data->meta_[0].length(skip_processing_meta_field.value().c_str())>0)
        {
            if (this->next()->putq(m1) == -1)
            {
                GERROR("CmrParametricMappingGadget::process, passing incoming image array on to next gadget");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        // -------------------------------------------------------------

        size_t encoding = (size_t)data->meta_[0].as_long("encoding", 0);
        GADGET_CHECK_RETURN(encoding<num_encoding_spaces_, GADGET_FAIL);

        std::string dataRole = std::string(data->meta_[0].as_str(GADGETRON_DATA_ROLE));

        size_t RO = data->data_.get_size(0);
        size_t E1 = data->data_.get_size(1);
        size_t E2 = data->data_.get_size(2);
        size_t CHA = data->data_.get_size(3);
        size_t N = data->data_.get_size(4);
        size_t S = data->data_.get_size(5);
        size_t SLC = data->data_.get_size(6);

        std::stringstream os;
        os << "_encoding_" << encoding << "_processing_" << process_called_times_;
        std::string str = os.str();

        // -------------------------------------------------------------

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array_complex(data->data_, debug_folder_full_path_ + "data" + str);
        }

        // -------------------------------------------------------------
        // if prep times are not read from protocol, find them from image header
        if (!imaging_prep_time_from_protocol.value())
        {
            this->prep_times_.resize(N, 0);

            size_t n;
            for (n = 0; n < N; n++)
            {
                this->prep_times_[n] = data->headers_(n).user_int[7] * 1e-3; // convert microsecond to ms
            }

            if (this->verbose.value())
            {
                GDEBUG_STREAM("Prep times are read from images ... ");
                for (n = 0; n < N; n++)
                {
                    GDEBUG_STREAM("Image " << n << " - " << this->prep_times_[n] << " ms ");
                }
            }
        }

        // -------------------------------------------------------------

        // calling the mapping
        IsmrmrdImageArray map, para, map_sd, para_sd;

        int status = this->perform_mapping(*data, map, para, map_sd, para_sd);

        if (status != GADGET_OK)
        {
            GWARN("CmrParametricMappingGadget::process, process incoming data failed ... ");

            // sending the incoming images
            if (this->next()->putq(m1) == -1)
            {
                GERROR("CmrParametricMappingGadget::process, passing incoming image array on to next gadget");
                return GADGET_FAIL;
            }
        }
        else
        {
            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(map.data_, debug_folder_full_path_ + "map" + str);
                gt_exporter_.export_array_complex(para.data_, debug_folder_full_path_ + "para" + str);
                gt_exporter_.export_array_complex(map_sd.data_, debug_folder_full_path_ + "map_sd" + str);
                gt_exporter_.export_array_complex(para_sd.data_, debug_folder_full_path_ + "para_sd" + str);
            }

            // sending the incoming images
            if (this->next()->putq(m1) == -1)
            {
                GERROR("CmrParametricMappingGadget::process, passing incoming image array on to next gadget");
                return GADGET_FAIL;
            }

            // ----------------------------------------------------------
            // fill in image header and meta
            // send out results
            // ----------------------------------------------------------
            if ( send_sd_map.value() && (this->fill_sd_header(map_sd) == GADGET_OK) )
            {
                Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>* cm2 = new Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>();
                *(cm2->getObjectPtr()) = map_sd;
                if (this->next()->putq(cm2) == -1)
                {
                    GERROR("CmrParametricMappingGadget::process, passing sd map on to next gadget");
                    return GADGET_FAIL;
                }
            }

            if (send_map.value() && (this->fill_map_header(map) == GADGET_OK))
            {
                Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>* cm1 = new Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>();
                *(cm1->getObjectPtr()) = map;
                if (this->next()->putq(cm1) == -1)
                {
                    GERROR("CmrParametricMappingGadget::process, passing map on to next gadget");
                    return GADGET_FAIL;
                }
            }
        }

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    int CmrParametricMappingGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "CmrParametricMappingGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
        }

        return GADGET_OK;
    }

    int CmrParametricMappingGadget::fill_map_header(IsmrmrdImageArray& map)
    {
        try
        {
            size_t RO  = map.data_.get_size(0);
            size_t E1  = map.data_.get_size(1);
            size_t E2  = map.data_.get_size(2);
            size_t CHA = map.data_.get_size(3);
            size_t N   = map.data_.get_size(4);
            size_t S   = map.data_.get_size(5);
            size_t SLC = map.data_.get_size(6);

            size_t e2, cha, n, s, slc;

            std::string lut = color_lut_map.value();
            if (this->field_strength_T_ > 2)
            {
                lut = color_lut_map_3T.value();
            }

            std::ostringstream ostr;
            ostr << "x" << (double)scaling_factor_map.value();
            std::string scalingStr = ostr.str();

            std::ostringstream ostr_unit;
            ostr_unit << std::setprecision(3) << 1.0f / scaling_factor_map.value() << "ms";
            std::string unitStr = ostr_unit.str();

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        size_t offset = n + s*N + slc*N*S;
                        map.meta_[offset].set(GADGETRON_IMAGE_SCALE_RATIO, (double)scaling_factor_map.value());
                        map.meta_[offset].set(GADGETRON_IMAGE_WINDOWCENTER, (long)(window_center_map.value()*scaling_factor_map.value()));
                        map.meta_[offset].set(GADGETRON_IMAGE_WINDOWWIDTH, (long)(window_width_map.value()*scaling_factor_map.value()));
                        map.meta_[offset].set(GADGETRON_IMAGE_COLORMAP, lut.c_str());

                        map.meta_[offset].set(GADGETRON_IMAGECOMMENT, map.meta_[offset].as_str(GADGETRON_DATA_ROLE));
                        map.meta_[offset].append(GADGETRON_IMAGECOMMENT, scalingStr.c_str());
                        map.meta_[offset].append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                    }
                }
            }

            Gadgetron::scal( (float)(scaling_factor_map.value()), map.data_);
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in CmrParametricMappingGadget::fill_map_header(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int CmrParametricMappingGadget::fill_sd_header(IsmrmrdImageArray& map_sd)
    {
        try
        {
            size_t RO = map_sd.data_.get_size(0);
            size_t E1 = map_sd.data_.get_size(1);
            size_t E2 = map_sd.data_.get_size(2);
            size_t CHA = map_sd.data_.get_size(3);
            size_t N = map_sd.data_.get_size(4);
            size_t S = map_sd.data_.get_size(5);
            size_t SLC = map_sd.data_.get_size(6);

            std::string lut = color_lut_sd_map.value();

            std::ostringstream ostr;
            ostr << "x" << (double)scaling_factor_sd_map.value();
            std::string scalingStr = ostr.str();

            std::ostringstream ostr_unit;
            ostr_unit << std::setprecision(3) << 1.0f / scaling_factor_sd_map.value() << "ms";
            std::string unitStr = ostr_unit.str();

            size_t e2, cha, n, s, slc;

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        size_t offset = n + s*N + slc*N*S;
                        map_sd.meta_[offset].set(GADGETRON_IMAGE_SCALE_RATIO, (double)scaling_factor_sd_map.value());
                        map_sd.meta_[offset].set(GADGETRON_IMAGE_WINDOWCENTER, (long)(window_center_sd_map.value()*scaling_factor_sd_map.value()));
                        map_sd.meta_[offset].set(GADGETRON_IMAGE_WINDOWWIDTH, (long)(window_width_sd_map.value()*scaling_factor_sd_map.value()));
                        map_sd.meta_[offset].set(GADGETRON_IMAGE_COLORMAP, lut.c_str());

                        map_sd.meta_[offset].set(GADGETRON_IMAGECOMMENT, map_sd.meta_[offset].as_str(GADGETRON_DATA_ROLE));
                        map_sd.meta_[offset].append(GADGETRON_IMAGECOMMENT, scalingStr.c_str());
                        map_sd.meta_[offset].append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                    }
                }
            }

            Gadgetron::scal((float)(scaling_factor_map.value()), map_sd.data_);
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in CmrParametricMappingGadget::fill_sd_header(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    bool CmrParametricMappingGadget::compute_mask_for_mapping(const hoNDArray<float>& mag, hoNDArray<float>& mask, float scale_factor)
    {
        try
        {
            size_t RO = mag.get_size(0);
            size_t E1 = mag.get_size(1);
            size_t SLC = mag.get_size(2);

            hoNDArray<float> mask_initial;
            mask_initial.create(RO, E1, SLC);
            Gadgetron::fill(mask_initial, (float)1);

            double std_thres = std_thres_masking.value();

            size_t n;
            for (n = 0; n < RO*E1*SLC; n++)
            {
                if (mag(n) < scale_factor*std_thres)
                {
                    mask_initial(n) = 0;
                }
            }

            mask = mask_initial;

            size_t obj_thres = 20;
            size_t bg_thres = 20;
            bool is_8_connected = true;

            size_t slc;
            for (slc = 0; slc < SLC; slc++)
            {
                hoNDArray<float> mask2D;
                mask2D.create(RO, E1, &mask_initial(0, 0, slc));

                hoNDArray<float> mask2D_clean;
                mask2D_clean.create(RO, E1, &mask(0, 0, slc));

                Gadgetron::bwlabel_clean_fore_and_background(mask2D, (float)1, (float)0, obj_thres, bg_thres, is_8_connected, mask2D_clean);
            }
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in CmrParametricMappingGadget::compute_mask_for_mapping(...) ... ");
            return false;
        }
    }

    // ----------------------------------------------------------------------------------------

    // GADGET_FACTORY_DECLARE(CmrParametricMappingGadget)

}
