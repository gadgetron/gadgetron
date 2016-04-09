
#include "FatWaterSeparationGadget.h"

namespace Gadgetron { 

FatWaterSeparationGadget::FatWaterSeparationGadget() : BaseClass()
{
}

FatWaterSeparationGadget::~FatWaterSeparationGadget()
{
}

int FatWaterSeparationGadget::process_config(ACE_Message_Block* mb)
{
    GADGET_CHECK_RETURN(BaseClass::process_config(mb)==GADGET_OK, GADGET_FAIL);

    ISMRMRD::IsmrmrdHeader h;
    try
    {
        deserialize(mb->rd_ptr(), h);
    }
    catch (...)
    {
        GDEBUG("Error parsing ISMRMRD Header");
        throw;
        return GADGET_FAIL;
    }

    GADGET_CHECK_RETURN(h.sequenceParameters.is_present(), GADGET_FAIL);
    GADGET_CHECK_RETURN(h.sequenceParameters->TE.is_present(), GADGET_FAIL);

    GADGET_CHECK_RETURN(h.sequenceParameters->TE->size() >= this->meas_max_idx_.contrast + 1, GADGET_FAIL);

    echo_times_.resize(this->meas_max_idx_.contrast + 1);

    size_t ind = 0;
    for (ind = 0; ind <= this->meas_max_idx_.contrast; ind++)
    {
        echo_times_[ind] = h.sequenceParameters->TE->operator[](ind);
    }

    if(this->verbose.value())
    {
        GDEBUG_STREAM("FatWaterSeparationGadget, find " << echo_times_.size() << " echoes ...");
        for (ind = 0; ind <= this->meas_max_idx_.contrast; ind++)
        {
            GDEBUG_STREAM("echo - " << ind << " - TE = " << echo_times_[ind] << "ms ");
        }
    }

    return GADGET_OK;
}

int FatWaterSeparationGadget::process_image_buffer(ImageBufferType& ori)
{
    GDEBUG_CONDITION_STREAM(verbose.value(), "FatWaterSeparationGadget::processImageBuffer(...) starts ... ");

    std::vector<std::string> processStr;
    std::vector<std::string> dataRole;

    boost::shared_ptr< std::vector<size_t> > dims = ori.get_dimensions();
    GDEBUG_CONDITION_STREAM(verbose.value(), "[Cha Slice E2 Con Phase Rep Set] = [" << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " << (*dims)[3] << " " << (*dims)[4]  << " " << (*dims)[5] << " " << (*dims)[6] << "]");

    // --------------------------------------------------------------------------------
    // ori
    // --------------------------------------------------------------------------------
    if (send_multi_echo_images.value())
    {
        GADGET_CHECK_RETURN(this->send_out_images(ori, processing_result_image_series_num.value(), processStr, dataRole)==GADGET_OK, GADGET_FAIL);
    }

    if (send_water_images.value() || send_fat_images.value() || send_t2_star_map.value() || send_field_map.value())
    {
        ImageBufferType water, fat, t2s_map, field_map;

        GADGET_CHECK_RETURN(this->perform_fat_water(ori, water, fat, t2s_map, field_map) == GADGET_OK, GADGET_FAIL);

        if (send_water_images.value())
        {
            processStr.clear();
            processStr.push_back(GADGETRON_IMAGE_WATER);

            dataRole.clear();
            dataRole.push_back(GADGETRON_IMAGE_WATER);

            GADGET_CHECK_RETURN(this->send_out_images(water, processing_result_image_series_num.value() +1, processStr, dataRole) == GADGET_OK, GADGET_FAIL);
        }

        if (send_fat_images.value())
        {
            processStr.clear();
            processStr.push_back(GADGETRON_IMAGE_FAT);

            dataRole.clear();
            dataRole.push_back(GADGETRON_IMAGE_FAT);

            GADGET_CHECK_RETURN(this->send_out_images(fat, processing_result_image_series_num.value() + 2, processStr, dataRole) == GADGET_OK, GADGET_FAIL);
        }

        if (send_t2_star_map.value())
        {
            processStr.clear();
            processStr.push_back(GADGETRON_IMAGE_T2STARMAP);

            dataRole.clear();
            dataRole.push_back(GADGETRON_IMAGE_T2STARMAP);

            GADGET_CHECK_RETURN(this->send_out_images(t2s_map, processing_result_image_series_num.value() + 3, processStr, dataRole) == GADGET_OK, GADGET_FAIL);
        }

        if (send_field_map.value())
        {
            processStr.clear();
            processStr.push_back(GADGETRON_IMAGE_FREQMAP);

            dataRole.clear();
            dataRole.push_back(GADGETRON_IMAGE_FREQMAP);

            GADGET_CHECK_RETURN(this->send_out_images(field_map, processing_result_image_series_num.value() + 4, processStr, dataRole) == GADGET_OK, GADGET_FAIL);
        }

        GADGET_CHECK_RETURN(this->release_image_buffer(water)     == GADGET_OK, GADGET_FAIL);
        GADGET_CHECK_RETURN(this->release_image_buffer(fat)       == GADGET_OK, GADGET_FAIL);
        GADGET_CHECK_RETURN(this->release_image_buffer(t2s_map)   == GADGET_OK, GADGET_FAIL);
        GADGET_CHECK_RETURN(this->release_image_buffer(field_map) == GADGET_OK, GADGET_FAIL);
    }

    GDEBUG_CONDITION_STREAM(verbose.value(), "FatWaterSeparationGadget::process(...) ends ... ");

    return GADGET_OK;
}

int FatWaterSeparationGadget::perform_fat_water(ImageBufferType& input, ImageBufferType& water, ImageBufferType& fat, ImageBufferType& t2s_map, ImageBufferType& field_map)
{
    size_t CHA = input.get_size(0);
    size_t SLC = input.get_size(1);
    size_t CON = input.get_size(2);
    size_t PHS = input.get_size(3);
    size_t REP = input.get_size(4);
    size_t SET = input.get_size(5);
    size_t AVE = input.get_size(6);

    std::vector<size_t> dim;
    input.get_dimensions(dim);

    dim[2] = 1; // only one contrast is outputted

    // set up outputs
    water.create(dim);
    fat.create(dim);
    t2s_map.create(dim);
    field_map.create(dim);

    // scaling and unit string for maps
    std::string scalingStr_T2StarMap, unitStr_T2StarMap;
    {
        std::ostringstream ostr;
        ostr << "x" << t2_star_scaling_factor.value();
        scalingStr_T2StarMap = ostr.str();
    }

    {
        std::ostringstream ostr;
        ostr << std::setprecision(3) << 1.0f / t2_star_scaling_factor.value() << "ms";
        unitStr_T2StarMap = ostr.str();
    }

    std::string scalingStr_FieldMap, unitStr_FieldMap, offsetStr_FieldMap;
    {
        std::ostringstream ostr;
        ostr << "x" << field_map_scaling_factor.value();
        scalingStr_FieldMap = ostr.str();
    }

    {
        std::ostringstream ostr;
        ostr << std::setprecision(3) << 1.0f / field_map_scaling_factor.value() << "Hz";
        unitStr_FieldMap = ostr.str();
    }

    {
        std::ostringstream ostr;
        ostr << "+" << field_map_offset.value();
        offsetStr_FieldMap = ostr.str();
    }

    if(verbose.value())
    {
        GDEBUG_STREAM("Scaling string for t2* map : " << scalingStr_T2StarMap);
        GDEBUG_STREAM("Unit string for t2* map : " << unitStr_T2StarMap);
        GDEBUG_STREAM("Scaling string for field map : " << scalingStr_FieldMap);
        GDEBUG_STREAM("Unit string for field map : " << unitStr_FieldMap);
        GDEBUG_STREAM("Offset string for field map : " << offsetStr_FieldMap);
    }

    size_t cha, slc, con, phs, rep, set, ave;

    ImageContainerType inputContainer, waterContainer, fatContainer, t2StarMapContainer, fieldMapContainer;

    std::vector<size_t> cols(SET, CON);
    inputContainer.create(cols, false);

    for ( cha=0; cha<CHA; cha++ )
    {
        for ( ave=0; ave<AVE; ave++ )
        {
            for ( slc=0; slc<SLC; slc++ )
            {
                for ( phs=0; phs<PHS; phs++ )
                {
                    for ( rep=0; rep<REP; rep++ )
                    {
                        for ( con=0; con<CON; con++ )
                        {
                            for ( set=0; set<SET; set++ )
                            {
                                inputContainer.set(input(cha, slc, con, phs, rep, set, ave), set, con);
                            }
                        }

                        /*if (!this->debug_folder_full_path_.empty())
                        {
                            std::ostringstream ostr;
                            ostr << "multi_echo_CHA" << cha << "_AVE" << ave << "_SLC" << slc << "_PHS" << phs << "_REP" << rep;
                            this->export_image_container(inputContainer, ostr.str());
                        }*/

                        GADGET_CHECK_RETURN(this->perform_fat_water(inputContainer, waterContainer, fatContainer, t2StarMapContainer, fieldMapContainer) == GADGET_OK, GADGET_FAIL);

                        /*if (!this->debug_folder_full_path_.empty())
                        {
                            std::ostringstream ostr;
                            ostr << "water_CHA" << cha << "_AVE" << ave << "_SLC" << slc << "_PHS" << phs << "_REP" << rep;
                            this->export_image_container(waterContainer, ostr.str());
                        }

                        if (!this->debug_folder_full_path_.empty())
                        {
                            std::ostringstream ostr;
                            ostr << "fat_CHA" << cha << "_AVE" << ave << "_SLC" << slc << "_PHS" << phs << "_REP" << rep;
                            this->export_image_container(fatContainer, ostr.str());
                        }

                        if (!this->debug_folder_full_path_.empty())
                        {
                            std::ostringstream ostr;
                            ostr << "t2star_map_CHA" << cha << "_AVE" << ave << "_SLC" << slc << "_PHS" << phs << "_REP" << rep;
                            this->export_image_container(t2StarMapContainer, ostr.str());
                        }

                        if (!this->debug_folder_full_path_.empty())
                        {
                            std::ostringstream ostr;
                            ostr << "field_map_CHA" << cha << "_AVE" << ave << "_SLC" << slc << "_PHS" << phs << "_REP" << rep;
                            this->export_image_container(fieldMapContainer, ostr.str());
                        }*/

                        for ( set=0; set<SET; set++ )
                        {
                            water(cha, slc, 0, phs, rep, set, ave) = &waterContainer(set, 0);
                            water(cha, slc, 0, phs, rep, set, ave)->header_ = input(cha, slc, 0, phs, rep, set, ave)->header_;
                            water(cha, slc, 0, phs, rep, set, ave)->attrib_ = input(cha, slc, 0, phs, rep, set, ave)->attrib_;
                            water(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_WATER);
                            water(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_WATER);
                            water(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_WATER);
                            water(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_WATER);

                            fat(cha, slc, 0, phs, rep, set, ave) = &fatContainer(set, 0);
                            fat(cha, slc, 0, phs, rep, set, ave)->header_ = input(cha, slc, 0, phs, rep, set, ave)->header_;
                            fat(cha, slc, 0, phs, rep, set, ave)->attrib_ = input(cha, slc, 0, phs, rep, set, ave)->attrib_;
                            fat(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_FAT);
                            fat(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_FAT);
                            fat(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_FAT);
                            fat(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_FAT);

                            t2s_map(cha, slc, 0, phs, rep, set, ave) = &t2StarMapContainer(set, 0);
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->header_ = input(cha, slc, 0, phs, rep, set, ave)->header_;
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_ = input(cha, slc, 0, phs, rep, set, ave)->attrib_;
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T2STARMAP);
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T2STARMAP);
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGE_SCALE_RATIO, t2_star_scaling_factor.value());
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGE_WINDOWCENTER, (long)(t2_star_window_center.value()*t2_star_scaling_factor.value()));
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGE_WINDOWWIDTH, (long)(t2_star_window_width.value()*t2_star_scaling_factor.value()));
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGE_COLORMAP, t2_star_color_map.value().c_str());
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T2STARMAP);

                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGECOMMENT, "GT");
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_T2STARMAP);
                            t2s_map(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGECOMMENT, unitStr_T2StarMap.c_str());

                            field_map(cha, slc, 0, phs, rep, set, ave) = &fieldMapContainer(set, 0);
                            field_map(cha, slc, 0, phs, rep, set, ave)->header_ = input(cha, slc, 0, phs, rep, set, ave)->header_;
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_ = input(cha, slc, 0, phs, rep, set, ave)->attrib_;
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_FREQMAP);
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_FREQMAP);
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGE_SCALE_RATIO, field_map_scaling_factor.value());
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGE_SCALE_OFFSET, field_map_offset.value());
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGE_WINDOWCENTER, (long)(field_map_window_center.value()*field_map_scaling_factor.value()));
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGE_WINDOWWIDTH, (long)(field_map_window_width.value()*field_map_scaling_factor.value()));
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGE_COLORMAP, field_map_color_map.value().c_str());
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_FREQMAP);

                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.set(GADGETRON_IMAGECOMMENT, "GT");
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_FREQMAP);
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGECOMMENT, unitStr_FieldMap.c_str());
                            field_map(cha, slc, 0, phs, rep, set, ave)->attrib_.append(GADGETRON_IMAGECOMMENT, offsetStr_FieldMap.c_str());
                        }

                        waterContainer.delete_data_on_destruct(false);
                        fatContainer.delete_data_on_destruct(false);
                        t2StarMapContainer.delete_data_on_destruct(false);
                        fieldMapContainer.delete_data_on_destruct(false);
                    }
                }
            }
        }
    }

    return GADGET_OK;
}

int FatWaterSeparationGadget::perform_fat_water(ImageContainerType& input, ImageContainerType& water, ImageContainerType& fat, ImageContainerType& t2s_map, ImageContainerType& field_map)
{
    try
    {
        std::vector<size_t> cols = input.cols();
        std::vector<size_t> cols_res(cols.size(), 1); // every row has one image as output

        water.create(cols_res, true);
        fat.create(cols_res, true);
        t2s_map.create(cols_res, true);
        field_map.create(cols_res, true);

        // for every row, call up fat water seperation on the multi-echo images
        size_t r, c;
        for (r = 0; r < cols.size(); r++)
        {
            water(r, 0) = input(r, 0);
            fat(r, 0) = input(r, 0);
            t2s_map(r, 0) = input(r, 0);
            field_map(r, 0) = input(r, 0);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in FatWaterSeparationGadget::perform_fat_water(ImageContainerType& input, ...) ... ");
        return GADGET_FAIL;
    }

    return GADGET_OK;
}

int FatWaterSeparationGadget::close(unsigned long flags)
{
    GDEBUG_CONDITION_STREAM(true, "FatWaterSeparationGadget - close(flags) : " << flags);

    if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(FatWaterSeparationGadget)

}
