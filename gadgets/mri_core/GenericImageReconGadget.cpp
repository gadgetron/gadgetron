
#include "GenericImageReconGadget.h"
#include <iomanip>

namespace Gadgetron {

    // ------------------------------------------------------------------------

    GenericImageReconGadgetBase::GenericImageReconGadgetBase()
    {
        image_series_num_ = 100;

        gt_timer_.set_timing_in_destruction(false);
        gt_timer_local_.set_timing_in_destruction(false);

        send_out_gfactor_map_ = false;
        send_out_snr_map_ = true;
        send_out_std_map_ = false;
        send_out_wrap_around_map_ = false;

        send_out_PSIR_as_real_ = false;
    }

    GenericImageReconGadgetBase::~GenericImageReconGadgetBase()
    {
    }

    hoMRImage<std::complex<float>, 3>* GenericImageReconGadgetBase::getImage3DFromImage2D(Image2DBufferType& ori, size_t cha, size_t slc, size_t con, size_t phs, size_t rep, size_t set, size_t ave)
    {
        Image2DType* pImage2D = &ori(cha, slc, 0, con, phs, rep, set, ave);
        GADGET_CHECK_THROW(pImage2D != NULL);

        size_t RO = pImage2D->get_size(0);
        size_t E1 = pImage2D->get_size(1);
        size_t E2 = ori.get_size(2);

        Image3DType* pImage3D = new Image3DType(RO, E1, E2);
        GADGET_CHECK_THROW(pImage3D != NULL);

        pImage3D->attrib_ = pImage2D->attrib_;

        size_t e2;
        for (e2 = 0; e2 < E2; e2++)
        {
            pImage2D = &ori(cha, slc, e2, con, phs, rep, set, ave);
            GADGET_CHECK_THROW(pImage2D != NULL);

            memcpy(pImage3D->begin() + e2 * RO*E1, pImage2D->begin(), sizeof(ValueType)*RO*E1);
        }

        return pImage3D;
    }

    bool GenericImageReconGadgetBase::getImage2DFromImage3D(Image3DType& image3D, Image2DBufferType& image2DBuf)
    {
        size_t RO = image3D.get_size(0);
        size_t E1 = image3D.get_size(1);
        size_t E2 = image3D.get_size(2);

        std::vector<size_t> dim(1);
        dim[0] = E2;
        image2DBuf.create(dim);

        size_t e2;
        for (e2 = 0; e2 < E2; e2++)
        {
            Image2DType* pImage2D = new Image2DType(RO, E1);
            GADGET_CHECK_RETURN_FALSE(pImage2D != NULL);

            memcpy(pImage2D->begin(), image3D.begin() + e2 * RO*E1, sizeof(ValueType)*RO*E1);

            image2DBuf(e2) = *pImage2D;
        }

        return true;
    }

    size_t GenericImageReconGadgetBase::computeSeriesImageNumber(ISMRMRD::ISMRMRD_ImageHeader& imheader, size_t nCHA, size_t cha, size_t nE2, size_t e2)
    {
        size_t nSET = meas_max_idx_.set + 1;
        size_t nREP = meas_max_idx_.repetition + 1;
        size_t nPHS = meas_max_idx_.phase + 1;
        size_t nSLC = meas_max_idx_.slice + 1;
        size_t nCON = meas_max_idx_.contrast + 1;
        if (nE2 == 0) nE2 = 1;

        size_t imageNum = imheader.average*nREP*nSET*nPHS*nCON*nSLC*nE2*nCHA
            + imheader.repetition*nSET*nPHS*nCON*nSLC*nE2*nCHA
            + imheader.set*nPHS*nCON*nSLC*nE2*nCHA
            + imheader.phase*nCON*nSLC*nE2*nCHA
            + imheader.contrast*nSLC*nE2*nCHA
            + imheader.slice*nE2*nCHA
            + e2 * nCHA
            + cha
            + 1;

        return imageNum;
    }

    void GenericImageReconGadgetBase::decorateImageHeader(ISMRMRD::ISMRMRD_ImageHeader& header, ISMRMRD::MetaContainer& attrib, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter, const std::vector<float>& windowWidth, bool resetImageCommentsParametricMaps, size_t slc, size_t SLC)
    {
        try
        {
            if (!dataRole.empty() && (dataRole[0] != GADGETRON_IMAGE_REGULAR))
            {
                std::string str;

                // data role
                bool isRealImage = false;
                bool isParametricMap = false;
                bool isParametricT1Map = false;
                bool isParametricT1SDMap = false;
                bool isParametricT2Map = false;
                bool isParametricT2SDMap = false;
                bool isParametricT2StarMap = false;
                bool isParametricR2StarMap = false;
                bool isParametricT2StarMaskMap = false;
                bool isParametricT2StarSDMap = false;
                bool isParametricT2StarAMap = false;
                bool isParametricT2StarTruncMap = false;
                bool isParametricPerfMap = false;
                bool isParametricPerfFlowMap = false;
                bool isParametricPerfInterstitialVolumeMap = false;
                bool isParametricPerfMTTMap = false;
                bool isParametricPerfBloodVolumeMap = false;
                bool isParametricPerfEMap = false;
                bool isParametricPerfPSMap = false;
                bool isParametricAIFGd = false;
                bool isParametricPerfGd = false;
                bool isParametricWater = false;
                bool isParametricFat = false;

                if (!dataRole.empty())
                {
                    size_t n;
                    for (n = 0; n < dataRole.size(); n++)
                    {
                        if (send_out_PSIR_as_real_ && dataRole[n] == GADGETRON_IMAGE_PSIR)
                        {
                            isRealImage = true;
                        }

                        if ((dataRole[n] == GADGETRON_IMAGE_T1MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_T1SDMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_T2MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_T2SDMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_T2STARMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_R2STARMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_T2STARMASKMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_T2STARSDMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_T2STARAMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_T2STARTRUNCMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_FREQMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_B1MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_FLIPANGLEMAP)
                            || (dataRole[n] == GADGETRON_IMAGE_PERF_MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_PERF_FLOW_MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_PERF_INTERSITITAL_VOLUME_MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_PERF_MEANTRANSITTIME_MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_PERF_VASCULAR_VOLUME_MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_PERF_Gd_Extraction_MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_PERF_PERMEABILITY_SURFACE_AREA_MAP)
                            || (dataRole[n] == GADGETRON_IMAGE_AIF_Gd_CONCENTRATION)
                            || (dataRole[n] == GADGETRON_IMAGE_PERF_Gd_CONCENTRATION)
                            || (dataRole[n] == GADGETRON_IMAGE_FAT)
                            || (dataRole[n] == GADGETRON_IMAGE_WATER)
                            )
                        {
                            isParametricMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_T1MAP)
                        {
                            isParametricT1Map = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_T1SDMAP)
                        {
                            isParametricT1SDMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_T2MAP)
                        {
                            isParametricT2Map = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_T2SDMAP)
                        {
                            isParametricT2SDMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_T2STARMAP)
                        {
                            isParametricT2StarMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_R2STARMAP)
                        {
                            isParametricR2StarMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_T2STARSDMAP)
                        {
                            isParametricT2StarSDMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_T2STARAMAP)
                        {
                            isParametricT2StarAMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_T2STARTRUNCMAP)
                        {
                            isParametricT2StarTruncMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_T2STARMASKMAP)
                        {
                            isParametricT2StarMaskMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_PERF_MAP)
                        {
                            isParametricPerfMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_PERF_FLOW_MAP)
                        {
                            isParametricPerfFlowMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_PERF_INTERSITITAL_VOLUME_MAP)
                        {
                            isParametricPerfInterstitialVolumeMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_PERF_MEANTRANSITTIME_MAP)
                        {
                            isParametricPerfMTTMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_PERF_VASCULAR_VOLUME_MAP)
                        {
                            isParametricPerfBloodVolumeMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_PERF_Gd_Extraction_MAP)
                        {
                            isParametricPerfEMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_PERF_PERMEABILITY_SURFACE_AREA_MAP)
                        {
                            isParametricPerfPSMap = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_AIF_Gd_CONCENTRATION)
                        {
                            isParametricAIFGd = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_PERF_Gd_CONCENTRATION)
                        {
                            isParametricPerfGd = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_FAT)
                        {
                            isParametricFat = true;
                        }

                        if (dataRole[n] == GADGETRON_IMAGE_WATER)
                        {
                            isParametricWater = true;
                        }
                    }

                    Gadgetron::append_ismrmrd_meta_values(attrib, GADGETRON_DATA_ROLE, dataRole);

                    std::vector<std::string> dataRoleAll;
                    Gadgetron::get_ismrmrd_meta_values(attrib, GADGETRON_DATA_ROLE, dataRoleAll);

                    // double check the image type
                    if (isRealImage)
                    {
                        header.image_type = ISMRMRD::ISMRMRD_IMTYPE_REAL;
                    }

                    // image comment
                    if (isParametricMap)
                    {
                        // reset the image comment for maps

                        std::vector<std::string> commentStr(dataRole.size() + 1);

                        commentStr[0] = "GT";
                        for (n = 0; n < dataRole.size(); n++)
                        {
                            commentStr[n + 1] = dataRole[n];
                        }

                        if (resetImageCommentsParametricMaps)
                        {
                            Gadgetron::set_ismrmrd_meta_values(attrib, GADGETRON_IMAGECOMMENT, commentStr);
                        }
                        else
                        {
                            for (n = 0; n < dataRole.size(); n++)
                            {
                                attrib.append(GADGETRON_IMAGECOMMENT, dataRole[n].c_str());
                            }
                        }

                        // get the scaling ratio
                        float scalingRatio = 1;
                        try
                        {
                            scalingRatio = (float)(attrib.as_double(GADGETRON_IMAGE_SCALE_RATIO, 0));

                            std::ostringstream ostr;
                            ostr << "x" << scalingRatio;
                            std::string scalingStr = ostr.str();
                            attrib.append(GADGETRON_IMAGECOMMENT, scalingStr.c_str());

                            if (isParametricT1Map || isParametricT1SDMap
                                || isParametricT2Map || isParametricT2SDMap
                                || isParametricT2StarMap || isParametricT2StarSDMap)
                            {
                                std::ostringstream ostr;
                                ostr << std::setprecision(3) << 1.0f / scalingRatio << "ms";
                                std::string unitStr = ostr.str();

                                attrib.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                            }

                            if (isParametricR2StarMap)
                            {
                                std::ostringstream ostr;
                                ostr << std::setprecision(3) << 1.0f / scalingRatio << "Hz";
                                std::string unitStr = ostr.str();

                                attrib.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                            }

                            if (isParametricPerfFlowMap || isParametricPerfMap || isParametricPerfPSMap)
                            {
                                std::ostringstream ostr;
                                ostr << std::setprecision(3) << 1.0f / scalingRatio << "ml/min/g";
                                std::string unitStr = ostr.str();

                                attrib.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                            }

                            if (isParametricPerfInterstitialVolumeMap)
                            {
                                std::ostringstream ostr;
                                ostr << std::setprecision(3) << 1.0f / scalingRatio << "ml/100g";
                                std::string unitStr = ostr.str();

                                attrib.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                            }

                            if (isParametricPerfMTTMap)
                            {
                                std::ostringstream ostr;
                                ostr << std::setprecision(3) << 1.0f / scalingRatio << "s";
                                std::string unitStr = ostr.str();

                                attrib.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                            }

                            if (isParametricPerfEMap)
                            {
                                std::ostringstream ostr;
                                ostr << std::setprecision(3) << 1.0f / scalingRatio;
                                std::string unitStr = ostr.str();

                                attrib.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                            }

                            if (isParametricPerfBloodVolumeMap)
                            {
                                std::ostringstream ostr;
                                ostr << std::setprecision(3) << 1.0f / scalingRatio << "ml/100g";
                                std::string unitStr = ostr.str();

                                attrib.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                            }

                            if (isParametricAIFGd || isParametricPerfGd)
                            {
                                std::ostringstream ostr;
                                ostr << std::setprecision(3) << 1.0f / scalingRatio << "MMOL/L";
                                std::string unitStr = ostr.str();

                                attrib.append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                            }
                        }
                        catch (...)
                        {
                            GWARN_STREAM("Image attrib does not have the scale ratio ...");
                            scalingRatio = 1;
                        }
                    }
                    else
                    {
                        bool hasSCC = false;
                        size_t n;
                        for (n = 0; n < processStr.size(); n++)
                        {
                            if (processStr[n] == GADGETRON_IMAGE_SURFACECOILCORRECTION)
                            {
                                hasSCC = true;
                            }
                        }

                        if (hasSCC)
                        {
                            double scalingRatio = (float)(attrib.as_double(GADGETRON_IMAGE_SCALE_RATIO, 0));

                            std::ostringstream ostr;
                            ostr << "x" << scalingRatio;
                            std::string scalingStr = ostr.str();

                            size_t numComment = attrib.length(GADGETRON_IMAGECOMMENT);
                            if (numComment == 0)
                            {
                                attrib.append(GADGETRON_IMAGECOMMENT, scalingStr.c_str());
                            }
                            else
                            {
                                size_t n;
                                for (n = 0; n < numComment; n++)
                                {
                                    std::string commentItem = attrib.as_str(GADGETRON_IMAGECOMMENT, n);
                                    if (commentItem.find("x") != std::string::npos)
                                    {
                                        const ISMRMRD::MetaValue& v = attrib.value(GADGETRON_IMAGECOMMENT, n);
                                        ISMRMRD::MetaValue& v2 = const_cast<ISMRMRD::MetaValue&>(v);
                                        v2 = scalingStr.c_str();
                                    }
                                }
                            }

                            if (attrib.length(GADGETRON_IMAGE_SCALE_OFFSET) > 0)
                            {
                                double offset = (float)(attrib.as_double(GADGETRON_IMAGE_SCALE_OFFSET, 0));

                                if (offset > 0)
                                {
                                    std::ostringstream ostr;
                                    ostr << "+" << offset;
                                    std::string offsetStr = ostr.str();

                                    attrib.append(GADGETRON_IMAGECOMMENT, offsetStr.c_str());
                                }
                            }
                        }

                        for (n = 0; n < dataRole.size(); n++)
                        {
                            attrib.append(GADGETRON_IMAGECOMMENT, dataRole[n].c_str());
                        }
                    }

                    // seq description
                    Gadgetron::append_ismrmrd_meta_values(attrib, GADGETRON_SEQUENCEDESCRIPTION, dataRole);

                    if(isParametricFat)
                    {
                        std::vector<std::string> dstr(1, GADGETRON_IMAGE_FAT);
                        Gadgetron::append_ismrmrd_meta_values(attrib, GADGETRON_SEQUENCEDESCRIPTION, dstr);
                    }

                    if (isParametricWater)
                    {
                        std::vector<std::string> dstr(1, GADGETRON_IMAGE_WATER);
                        Gadgetron::append_ismrmrd_meta_values(attrib, GADGETRON_SEQUENCEDESCRIPTION, dstr);
                    }
                }

                // image processing history
                if (!processStr.empty())
                {
                    size_t n;
                    for (n = 0; n < processStr.size(); n++)
                    {
                        attrib.append(GADGETRON_IMAGEPROCESSINGHISTORY, processStr[n].c_str());
                    }
                }

                if (windowCenter.size() == SLC && windowWidth.size() == SLC)
                {
                    if(windowCenter[slc]>0) attrib.set(GADGETRON_IMAGE_WINDOWCENTER, (long)windowCenter[slc]);
                    if(windowWidth[slc]>0) attrib.set(GADGETRON_IMAGE_WINDOWWIDTH, (long)windowWidth[slc]);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Exceptions happened in decorateImageHeader(...) ... ");
        }
    }

    bool GenericImageReconGadgetBase::convertContainer3Dto2D(const ImageContainer3DType& input, std::vector<size_t>& e2, ImageContainer2DType& output)
    {
        try
        {
            std::vector<size_t> cols = input.cols();
            std::vector<size_t> cols_output = cols;

            e2.resize(cols.size(), 0);

            size_t r, c;
            for (r = 0; r < cols.size(); r++)
            {
                GADGET_CHECK_RETURN_FALSE(input.has_identical_dimensions(r));
            }

            for (r = 0; r < cols.size(); r++)
            {
                size_t num2D = 0;
                for (c = 0; c < cols[r]; c++)
                {
                    size_t E2 = input(r, c).get_size(2);
                    e2[r] = E2;
                    num2D += E2;
                }

                cols_output[r] = num2D;
            }

            output.create(cols_output, true);

            std::vector<size_t> dim2D(2);

            for (r = 0; r < cols.size(); r++)
            {
                size_t num2D = 0;
                for (c = 0; c < cols[r]; c++)
                {
                    size_t RO = input(r, c).get_size(0);
                    size_t E1 = input(r, c).get_size(1);
                    size_t E2 = input(r, c).get_size(2);

                    dim2D[0] = RO;
                    dim2D[1] = E1;

                    for (size_t e2in3D = 0; e2in3D < E2; e2in3D++)
                    {
                        output(r, num2D + e2in3D).create(dim2D);
                        memcpy(output(r, num2D + e2in3D).begin(), &(input(r, c)(0, 0, e2in3D)), sizeof(T)*RO*E1);

                        output(r, num2D + e2in3D).header_ = input(r, c).header_;
                        output(r, num2D + e2in3D).attrib_ = input(r, c).attrib_;

                        output(r, num2D + e2in3D).header_.matrix_size[2] = 1;
                    }

                    num2D += E2;
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadgetBase::convertContainer3Dto2D(...) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadgetBase::convertContainer2Dto3D(const ImageContainer2DType& input, const std::vector<size_t>& e2, ImageContainer3DType& output)
    {
        try
        {
            std::vector<size_t> cols = input.cols();
            std::vector<size_t> cols_output(cols.size(), 0);

            GADGET_CHECK_RETURN_FALSE(e2.size() == cols.size());

            size_t r, c, outputc;
            for (r = 0; r < cols.size(); r++)
            {
                size_t num3D = cols[r] / e2[r];
                cols_output[r] = num3D;
            }

            output.create(cols_output, true);

            std::vector<size_t> dim3D(3);

            for (r = 0; r < cols.size(); r++)
            {
                size_t RO = input(r, 0).get_size(0);
                size_t E1 = input(r, 0).get_size(1);

                dim3D[0] = RO;
                dim3D[1] = E1;
                dim3D[2] = e2[r];

                size_t ind2D = 0;
                for (outputc = 0; outputc < cols_output[r]; outputc++)
                {
                    output(r, outputc).create(dim3D);
                    output(r, outputc).header_ = input(r, 0).header_;
                    output(r, outputc).attrib_ = input(r, 0).attrib_;

                    output(r, outputc).header_.matrix_size[2] = dim3D[2];

                    for (c = ind2D; c < ind2D + e2[r]; c++)
                    {
                        memcpy(output(r, outputc).begin() + (c - ind2D)*RO*E1, input(r, ind2D).begin(), sizeof(T)*RO*E1);
                    }

                    ind2D += e2[r];
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadgetBase::convertContainer2Dto3D(...) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer2D(ImageContainer2DType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<ValueType> outArray;

            for (r = 0; r < R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(outArray, debug_folder_full_path_ + ostr.str());
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer2D(ImageContainer2DMagType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<T> outArray;

            for (r = 0; r < R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array(outArray, debug_folder_full_path_ + ostr.str());
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer2D(ImageContainer2DDoubleType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray< std::complex<double> > outArray;

            for (r = 0; r < R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(outArray, debug_folder_full_path_ + ostr.str());
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer2D(ImageContainer2DMagDoubleType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<double> outArray;

            for (r = 0; r < R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array(outArray, debug_folder_full_path_ + ostr.str());
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer3D(ImageContainer3DType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r, c;
            for (r = 0; r < R; r++)
            {
                for (c = 0; c < input.cols(r); c++)
                {
                    std::ostringstream ostr;
                    ostr << prefix << "_" << r << "_" << c;

                    if (!debug_folder_full_path_.empty()) gt_exporter_.export_image_complex(input(r, c), debug_folder_full_path_ + ostr.str());
                }
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer3D(ImageContainer3DMagType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r, c;
            for (r = 0; r < R; r++)
            {
                for (c = 0; c < input.cols(r); c++)
                {
                    std::ostringstream ostr;
                    ostr << prefix << "_" << r << "_" << c;

                    if (!debug_folder_full_path_.empty()) gt_exporter_.export_image(input(r, c), debug_folder_full_path_ + ostr.str());
                }
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer2D(NDImageContainer2DType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<ValueType> outArray;

            for (r = 0; r < R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(outArray, debug_folder_full_path_ + ostr.str());
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer2D(NDImageContainer2DMagType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<T> outArray;

            for (r = 0; r < R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array(outArray, debug_folder_full_path_ + ostr.str());
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer2D(NDImageContainer2DDoubleType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray< std::complex<double> > outArray;

            for (r = 0; r < R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(outArray, debug_folder_full_path_ + ostr.str());
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer2D(NDImageContainer2DMagDoubleType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<double> outArray;

            for (r = 0; r < R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array(outArray, debug_folder_full_path_ + ostr.str());
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer3D(NDImageContainer3DType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r, c;
            for (r = 0; r < R; r++)
            {
                for (c = 0; c < input.cols(r); c++)
                {
                    std::ostringstream ostr;
                    ostr << prefix << "_" << r << "_" << c;

                    if (!debug_folder_full_path_.empty()) gt_exporter_.export_image_complex(input(r, c), debug_folder_full_path_ + ostr.str());
                }
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::exportImageContainer3D(NDImageContainer3DMagType& input, const std::string& prefix)
    {
        if (!debug_folder_full_path_.empty())
        {
            size_t R = input.rows();

            size_t r, c;
            for (r = 0; r < R; r++)
            {
                for (c = 0; c < input.cols(r); c++)
                {
                    std::ostringstream ostr;
                    ostr << prefix << "_" << r << "_" << c;

                    if (!debug_folder_full_path_.empty()) gt_exporter_.export_image(input(r, c), debug_folder_full_path_ + ostr.str());
                }
            }
        }

        return true;
    }

    bool GenericImageReconGadgetBase::storeGFactorMap(const Image2DType& gmap)
    {
        try
        {
            if (gfactor_buf_.empty())
            {
                gfactor_buf_.push_back(gmap);
                return true;
            }

            size_t N = gfactor_buf_.size();

            size_t n;
            for (n = 0; n < N; n++)
            {
                long nSLC = gfactor_buf_[n].header_.slice;
                long slc = gmap.header_.slice;

                if (nSLC == slc)
                {
                    return true;
                }
            }

            gfactor_buf_.push_back(gmap);
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GenericImageReconGadget::storeGFactorMap(const Image2DType& gmap) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadgetBase::findGFactorMap(size_t slc, Image2DType& gmap)
    {
        if (gfactor_buf_.empty()) return false;

        size_t N = gfactor_buf_.size();

        size_t n;
        for (n = 0; n < N; n++)
        {
            long nSLC = gfactor_buf_[n].header_.slice;

            if (nSLC == slc)
            {
                gmap = gfactor_buf_[n];
                return true;
            }
        }

        return false;
    }

    bool GenericImageReconGadgetBase::storeGFactorMap(const Image3DType& gmap)
    {
        try
        {
            if (gfactor_buf_3D_.empty())
            {
                gfactor_buf_3D_.push_back(gmap);
                return true;
            }

            size_t N = gfactor_buf_3D_.size();

            size_t n;
            for (n = 0; n < N; n++)
            {
                long nSLC = gfactor_buf_3D_[n].header_.slice;
                long slc = gmap.header_.slice;

                if (nSLC == slc)
                {
                    return true;
                }
            }

            gfactor_buf_3D_.push_back(gmap);
        }
        catch (...)
        {
            GERROR_STREAM("Error happened in GenericImageReconGadget::storeGFactorMap(const Image3DType& gmap) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadgetBase::findGFactorMap(size_t slc, Image3DType& gmap)
    {
        if (gfactor_buf_3D_.empty()) return false;

        size_t N = gfactor_buf_3D_.size();

        size_t n;
        for (n = 0; n < N; n++)
        {
            long nSLC = gfactor_buf_3D_[n].header_.slice;

            if (nSLC == slc)
            {
                gmap = gfactor_buf_3D_[n];
                return true;
            }
        }

        return false;
    }

    bool GenericImageReconGadgetBase::set_pixel_value_range(hoNDArray<ValueType>& im, T min_V, T max_V, bool ignore_empty_image)
    {
        size_t N = im.get_number_of_elements();

        if (ignore_empty_image)
        {
            T v = Gadgetron::nrm2(im);
            if (v < FLT_EPSILON) return true;
        }

        bool has_max = false;
        bool has_min = false;

        size_t n;
        for (n = 0; n < N; n++)
        {
            if (im(n).real() > max_V)
            {
                im(n) = max_V;
                has_max = true;
            }

            if (im(n).real() < min_V)
            {
                im(n) = min_V;
                has_min = true;
            }
        }

        if (!has_max) im(0) = max_V;

        return true;
    }

    bool GenericImageReconGadgetBase::set_pixel_value_range(hoNDArray<T>& im, T min_V, T max_V, bool ignore_empty_image)
    {
        size_t N = im.get_number_of_elements();

        if (ignore_empty_image)
        {
            T v = Gadgetron::nrm2(im);
            if (v < FLT_EPSILON) return true;
        }

        bool has_max = false;
        bool has_min = false;

        size_t n;
        for (n = 0; n < N; n++)
        {
            if (im(n) > max_V)
            {
                im(n) = max_V;
                has_max = true;
            }

            if (im(n) < min_V)
            {
                im(n) = min_V;
                has_min = true;
            }
        }

        if (!has_max) im(0) = max_V;

        return true;
    }

    // ------------------------------------------------------------------------

    GenericImageReconGadget::GenericImageReconGadget() : GenericImageReconGadgetBase()
    {

    }

    GenericImageReconGadget::~GenericImageReconGadget()
    {

    }

    bool GenericImageReconGadget::readParameters()
    {
        try
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "------> GenericImageReconGadget parameters <------");

            image_series_num_ = image_series_num.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "image_series_num_ is " << image_series_num_);

            GDEBUG_CONDITION_STREAM(verbose.value(), "-----------------------------------------------");

            send_out_gfactor_map_ = send_out_gfactor_map.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_out_gfactor_map_ is " << send_out_gfactor_map_);

            send_out_snr_map_ = send_out_snr_map.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_out_snr_map_ is " << send_out_snr_map_);

            send_out_std_map_ = send_out_std_map.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_out_std_map_ is " << send_out_std_map_);

            send_out_wrap_around_map_ = send_out_wrap_around_map.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_out_wrap_around_map_ is " << send_out_wrap_around_map_);

            send_out_PSIR_as_real_ = send_out_PSIR_as_real.value();
            GDEBUG_CONDITION_STREAM(verbose.value(), "send_out_PSIR_as_real_ is " << send_out_PSIR_as_real_);

            GDEBUG_CONDITION_STREAM(verbose.value(), "-----------------------------------------------");
        }
        catch (...)
        {
            GERROR_STREAM("Errors in GenericImageReconGadget::readParameters() ... ");
            return false;
        }

        return true;
    }

    int GenericImageReconGadget::process_config(ACE_Message_Block* mb)
    {
        // read in parameters from the xml
        GADGET_CHECK_RETURN(this->readParameters(), GADGET_FAIL);

        if (!debug_folder.value().empty())
        {
            Gadgetron::get_debug_folder_path(debug_folder.value(), debug_folder_full_path_);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is " << debug_folder_full_path_);

            // Create debug folder if necessary
            boost::filesystem::path boost_folder_path(debug_folder_full_path_);
            try
            {
                boost::filesystem::create_directories(boost_folder_path);
            }
            catch (...)
            {
                GERROR("Error creating the debug folder.\n");
                return false;
            }
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is not set ... ");
        }

        ISMRMRD::IsmrmrdHeader h;
        try {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...) {
            GDEBUG("Error parsing ISMRMRD Header");
            throw;
            return GADGET_FAIL;
        }

        // seq object
        if (h.encoding.size() != 1)
        {
            GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
        }

        GADGET_CATCH_THROW(find_encoding_limits(h, meas_max_idx_, verbose.value()));

        // ---------------------------------------------------------
        // check interpolation is on or off
        // ---------------------------------------------------------

        find_matrix_size_encoding(h, matrix_size_encoding_);
        find_matrix_size_recon(h, matrix_size_recon_);

        GDEBUG_CONDITION_STREAM(verbose.value(), "matrix size for encoding : [" << matrix_size_encoding_[0] << " " << matrix_size_encoding_[1] << " " << matrix_size_encoding_[2] << "]");
        GDEBUG_CONDITION_STREAM(verbose.value(), "matrix size for recon : [" << matrix_size_recon_[0] << " " << matrix_size_recon_[1] << " " << matrix_size_recon_[2] << "]");

        if (matrix_size_recon_[0] > 0.75*matrix_size_encoding_[0])
        {
            inteprolation_on_ = true;
            GDEBUG_CONDITION_STREAM(verbose.value(), "inteprolation_on_ : " << inteprolation_on_);
        }

        return GADGET_OK;
    }

    int GenericImageReconGadget::process(ACE_Message_Block *mb)
    {
        GadgetContainerMessage< hoNDObjectArray< hoMRImage<ValueType, 2> > >* m1 = nullptr;
        GadgetContainerMessage< hoNDObjectArray< hoMRImage<ValueType, 3> > >* m2 = nullptr;
        GadgetContainerMessage< ISMRMRD::ImageHeader >* m3 = nullptr;

        if (m1 = AsContainerMessage<Image2DBufferType>(mb))
        {
            return this->process2D(m1);
        }
        else if (m2 = AsContainerMessage<Image3DBufferType>(mb))
        {
            return this->process3D(m2);
        }
        else if (m3 = AsContainerMessage<ISMRMRD::ImageHeader>(mb))
        {
            return this->process_image(m3);
        }
        else if (pass_on_undesired_data_)
        {
            return this->next()->putq(mb);
        }
        else
        {
            GERROR("GenericImageReconGadget Gadget received wrong input data type");
            return -1;
        }
    }

    int GenericImageReconGadget::process_image(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1)
    {
        return this->next()->putq(m1);
    }

    int GenericImageReconGadget::process2D(GadgetContainerMessage<Image2DBufferType>* m1)
    {
        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        Image2DBufferType ori = *m1->getObjectPtr();
        m1->release();

        if (ori.get_number_of_elements() == 1)
        {
            size_t num = ori(0).attrib_.length(GADGETRON_DATA_ROLE);
            GADGET_CHECK_RETURN(num > 0, GADGET_FAIL);

            dataRole.resize(num);

            for (size_t ii = 0; ii < num; ii++)
            {
                dataRole[ii] = std::string(ori(0).attrib_.as_str(GADGETRON_DATA_ROLE, ii));
            }

            if (dataRole[0] == GADGETRON_IMAGE_GFACTOR)
            {
                GADGET_CHECK_RETURN(this->storeGFactorMap(ori(0)), GADGET_FAIL);

                if (send_out_gfactor_map_)
                {
                    GDEBUG_STREAM("Image Recon, send out received gfactor map ...");

                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 200, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 200, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                }

                return GADGET_OK;
            }
            else if (dataRole[0] == GADGETRON_IMAGE_SNR_MAP)
            {
                if (send_out_snr_map_)
                {
                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 201, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 201, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                }

                return GADGET_OK;
            }
            else if (dataRole[0] == GADGETRON_IMAGE_STD_MAP)
            {
                if (send_out_std_map_)
                {
                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 202, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 202, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                }

                return GADGET_OK;
            }
            else if (dataRole[0] == GADGETRON_IMAGE_WRAPAROUNDMAP)
            {
                if (send_out_wrap_around_map_)
                {
                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 203, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 203, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                }

                return GADGET_OK;
            }

            num = ori(0).attrib_.length(GADGETRON_PASS_IMMEDIATE);
            if (num > 0)
            {
                long pass_image_immediately = ori(0).attrib_.as_long(GADGETRON_PASS_IMMEDIATE, 0);
                if (pass_image_immediately)
                {
                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 1004, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 1004, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                    return GADGET_OK;
                }
            }
        }

        this->processImageBuffer(ori);

        this->releaseImageBuffer(ori);

        return GADGET_OK;
    }

    int GenericImageReconGadget::process3D(GadgetContainerMessage<Image3DBufferType>* m1)
    {
        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        Image3DBufferType ori = *m1->getObjectPtr();
        m1->release();

        if (ori.get_number_of_elements() == 1)
        {
            size_t num = ori(0).attrib_.length(GADGETRON_DATA_ROLE);
            GADGET_CHECK_RETURN(num > 0, GADGET_FAIL);

            dataRole.resize(num);

            for (size_t ii = 0; ii < num; ii++)
            {
                dataRole[ii] = std::string(ori(0).attrib_.as_str(GADGETRON_DATA_ROLE, ii));
            }

            if (dataRole[0] == GADGETRON_IMAGE_GFACTOR)
            {
                GADGET_CHECK_RETURN(this->storeGFactorMap(ori(0)), GADGET_FAIL);

                if (send_out_gfactor_map_)
                {
                    GDEBUG_STREAM("Image Recon, send out received gfactor map ...");

                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 200, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 200, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                }

                return GADGET_OK;
            }
            else if (dataRole[0] == GADGETRON_IMAGE_SNR_MAP)
            {
                if (send_out_snr_map_)
                {
                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 201, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 201, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                }

                return GADGET_OK;
            }
            else if (dataRole[0] == GADGETRON_IMAGE_STD_MAP)
            {
                if (send_out_std_map_)
                {
                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 202, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 202, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                }

                return GADGET_OK;
            }
            else if (dataRole[0] == GADGETRON_IMAGE_WRAPAROUNDMAP)
            {
                if (send_out_wrap_around_map_)
                {
                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 203, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 203, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                }

                return GADGET_OK;
            }

            num = ori(0).attrib_.length(GADGETRON_PASS_IMMEDIATE);
            if (num > 0)
            {
                long pass_image_immediately = ori(0).attrib_.as_long(GADGETRON_PASS_IMMEDIATE, 0);
                if (pass_image_immediately)
                {
                    dataRole.clear();
                    if (send_out_image_array.value())
                    {
                        GADGET_CHECK_RETURN(this->sendOutImageBuffer(ori, image_series_num_ + 1004, processStr, dataRole), GADGET_FAIL);
                    }
                    else
                    {
                        GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_ + 1004, processStr, dataRole), GADGET_FAIL);
                    }
                    GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                    return GADGET_OK;
                }
            }
        }

        this->processImageBuffer(ori);

        this->releaseImageBuffer(ori);

        return GADGET_OK;
    }

    int GenericImageReconGadget::processImageBuffer(Image2DBufferType& ori)
    {
        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        size_t RO = ori(0).get_size(0);
        size_t E1 = ori(0).get_size(1);
        size_t E2 = ori(0).get_size(2);

        std::vector<size_t> dims = ori.get_dimensions();
        GDEBUG_CONDITION_STREAM(verbose.value(), "[RO E1 E2 Cha Slice Con Phase Rep Set Ave] = ["
            << RO << " " << E1 << " " << E2 << " "
            << dims[0] << " " << dims[1] << " " << dims[2] << " "
            << dims[3] << " " << dims[4] << " " << dims[5] << " "
            << dims[6] << " " << dims[7] << "]");

        this->sendOutImages(ori, image_series_num_++, processStr, dataRole);

        return GADGET_OK;
    }

    int GenericImageReconGadget::processImageBuffer(Image3DBufferType& ori)
    {
        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        size_t RO = ori(0).get_size(0);
        size_t E1 = ori(0).get_size(1);
        size_t E2 = ori(0).get_size(2);

        std::vector<size_t> dims = ori.get_dimensions();
        GDEBUG_CONDITION_STREAM(verbose.value(), "[RO E1 E2 Cha Slice Con Phase Rep Set Ave] = ["
            << RO << " " << E1 << " " << E2 << " "
            << dims[0] << " " << dims[1] << " " << dims[2] << " "
            << dims[3] << " " << dims[4] << " " << dims[5] << " "
            << dims[6] << " " << dims[7] << "]");

        this->sendOutImages(ori, image_series_num_++, processStr, dataRole);

        return GADGET_OK;
    }

    int GenericImageReconGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GenericImageReconGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
            std::string procTime;
            Gadgetron::get_current_moment(procTime);

            GDEBUG_STREAM("* ============================================================================== *");
            GDEBUG_STREAM("---> Image recon phase, Current processing time : " << procTime << " <---");
            GDEBUG_STREAM("* ============================================================================== *");
        }

        return GADGET_OK;
    }

    bool GenericImageReconGadget::sendOutImages(Image2DBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole,
        const std::vector<float>& windowCenter, const std::vector<float>& windowWidth, bool resetImageCommentsParametricMaps, Gadget* anchor)
    {
        try
        {
            size_t CHA = images.get_size(0);
            size_t SLC = images.get_size(1);
            size_t CON = images.get_size(2);
            size_t PHS = images.get_size(3);
            size_t REP = images.get_size(4);
            size_t SET = images.get_size(5);
            size_t AVE = images.get_size(6);

            std::string dataRoleString;
            if (!dataRole.empty())
            {
                std::ostringstream ostr;
                for (size_t n = 0; n < dataRole.size(); n++)
                {
                    ostr << dataRole[n] << " - ";
                }

                dataRoleString = ostr.str();
            }

            int series_no = seriesNum;
            if (seriesNum < 0)
            {
                for (size_t ii = 0; ii < images.get_number_of_elements(); ii++)
                {
                    if (images(ii).header_.image_series_index >= 0)
                    {
                        series_no = images(ii).header_.image_series_index;
                        break;
                    }
                }
            }

            GDEBUG_CONDITION_STREAM(verbose.value(), "--> GenericImageReconGadget, sending out " << dataRoleString << " images for series " << series_no << ", array boundary [CHA SLC CON PHS REP SET AVE] = ["
                << CHA << " " << SLC << " " << CON << " " << PHS << " " << REP << " " << SET << " " << AVE << "] ");

            size_t ave(0), set(0), rep(0), phs(0), con(0), slc(0), cha(0);
            std::vector<size_t> dim3D(3);

            for (ave = 0; ave < AVE; ave++)
            {
                for (set = 0; set < SET; set++)
                {
                    for (rep = 0; rep < REP; rep++)
                    {
                        for (phs = 0; phs < PHS; phs++)
                        {
                            for (con = 0; con < CON; con++)
                            {
                                for (slc = 0; slc < SLC; slc++)
                                {
                                    for (cha = 0; cha < CHA; cha++)
                                    {
                                        hoMRImage<ValueType, 2>* pImage = &images(cha, slc, con, phs, rep, set, ave);
                                        if (pImage != NULL && pImage->get_number_of_elements()>0)
                                        {
                                            T v = Gadgetron::nrm2(*pImage);
                                            if (v < FLT_EPSILON) continue; // do not send out empty image

                                            Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>();
                                            Gadgetron::GadgetContainerMessage<ImgArrayType>* cm2 = new Gadgetron::GadgetContainerMessage<ImgArrayType>();
                                            Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>* cm3 = new Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>();

                                            try
                                            {
                                                cm1->cont(cm2);
                                                cm2->cont(cm3);

                                                // set the ISMRMRD image header
                                                memcpy(cm1->getObjectPtr(), &pImage->header_, sizeof(ISMRMRD::ISMRMRD_ImageHeader));

                                                long long imageNum(0);
                                                if (pImage->attrib_.length(GADGETRON_IMAGENUMBER) == 0)
                                                {
                                                    imageNum = this->computeSeriesImageNumber(*cm1->getObjectPtr(), CHA, cha, 1, 0);
                                                    cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                    pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                }
                                                else
                                                {
                                                    imageNum = pImage->attrib_.as_long(GADGETRON_IMAGENUMBER);
                                                    if (imageNum > 0)
                                                    {
                                                        cm1->getObjectPtr()->image_index = imageNum;
                                                    }
                                                    else
                                                    {
                                                        imageNum = this->computeSeriesImageNumber(*cm1->getObjectPtr(), CHA, cha, 1, 0);
                                                        cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                        pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                    }
                                                }

                                                cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;

                                                if(add_original_series_num.value())
                                                {
                                                    if (seriesNum >= 0) cm1->getObjectPtr()->image_series_index += seriesNum;
                                                }
                                                else
                                                {
                                                    if (seriesNum >= 0) cm1->getObjectPtr()->image_series_index = seriesNum;
                                                }

                                                // set the image data
                                                size_t RO = pImage->get_size(0);
                                                size_t E1 = pImage->get_size(1);
                                                size_t E2 = pImage->get_size(2);

                                                dim3D[0] = RO;
                                                dim3D[1] = E1;
                                                dim3D[2] = E2;

                                                cm1->getObjectPtr()->matrix_size[0] = RO;
                                                cm1->getObjectPtr()->matrix_size[1] = E1;
                                                cm1->getObjectPtr()->matrix_size[2] = E2;

                                                cm2->getObjectPtr()->create(dim3D);
                                                memcpy(cm2->getObjectPtr()->get_data_ptr(), pImage->get_data_ptr(), pImage->get_number_of_bytes());

                                                // set the attributes
                                                *cm3->getObjectPtr() = pImage->attrib_;

                                                this->decorateImageHeader(*cm1->getObjectPtr(), *cm3->getObjectPtr(), seriesNum, processStr, dataRole, windowCenter, windowWidth, resetImageCommentsParametricMaps, slc, SLC);

                                                if (anchor != NULL)
                                                {
                                                    if (anchor->putq(cm1) < 0)
                                                    {
                                                        cm1->release();
                                                        return false;
                                                    }
                                                }
                                                else
                                                {
                                                    if (this->next()->putq(cm1) < 0)
                                                    {
                                                        cm1->release();
                                                        return false;
                                                    }
                                                }
                                            }
                                            catch (...)
                                            {
                                                cm1->release();
                                                throw;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadget::sendOutImages(images, seriesNum, processStr, dataRole) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadget::sendOutImageBuffer(Image2DBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole,
        const std::vector<float>& windowCenter, const std::vector<float>& windowWidth, bool resetImageCommentsParametricMaps, Gadget* anchor)
    {
        try
        {
            size_t CHA = images.get_size(0);
            size_t SLC = images.get_size(1);
            size_t CON = images.get_size(2);
            size_t PHS = images.get_size(3);
            size_t REP = images.get_size(4);
            size_t SET = images.get_size(5);
            size_t AVE = images.get_size(6);

            std::string dataRoleString;
            if (!dataRole.empty())
            {
                std::ostringstream ostr;
                for (size_t n = 0; n < dataRole.size(); n++)
                {
                    ostr << dataRole[n] << " - ";
                }

                dataRoleString = ostr.str();
            }

            GDEBUG_CONDITION_STREAM(verbose.value(), "--> GenericImageReconGadget, sending out image array " << dataRoleString << " images for series " << seriesNum << ", array boundary [CHA SLC CON PHS REP SET AVE] = ["
                << CHA << " " << SLC << " " << CON << " " << PHS << " " << REP << " " << SET << " " << AVE << "] ");

            size_t ave(0), set(0), rep(0), phs(0), con(0), slc(0), cha(0);
            std::vector<size_t> dim3D(3);

            for (ave = 0; ave < AVE; ave++)
            {
                for (set = 0; set < SET; set++)
                {
                    for (rep = 0; rep < REP; rep++)
                    {
                        for (phs = 0; phs < PHS; phs++)
                        {
                            for (con = 0; con < CON; con++)
                            {
                                for (slc = 0; slc < SLC; slc++)
                                {
                                    for (cha = 0; cha < CHA; cha++)
                                    {
                                        hoMRImage<ValueType, 2>* pImage = &images(cha, slc, con, phs, rep, set, ave);
                                        if (pImage != NULL && pImage->get_number_of_elements() > 0)
                                        {
                                            T v = Gadgetron::nrm2(*pImage);
                                            if (v < FLT_EPSILON) continue; // do not send out empty image

                                            try
                                            {
                                                long long imageNum(0);
                                                if (pImage->attrib_.length(GADGETRON_IMAGENUMBER) == 0)
                                                {
                                                    imageNum = this->computeSeriesImageNumber(pImage->header_, CHA, cha, 1, 0);
                                                    pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                }
                                                else
                                                {
                                                    imageNum = pImage->attrib_.as_long(GADGETRON_IMAGENUMBER);
                                                    if (imageNum <= 0)
                                                    {
                                                        imageNum = this->computeSeriesImageNumber(pImage->header_, CHA, cha, 1, 0);
                                                        pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                    }
                                                }

                                                pImage->header_.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                                                if (add_original_series_num.value())
                                                {
                                                    if (seriesNum >= 0) pImage->header_.image_series_index += seriesNum;
                                                }
                                                else
                                                {
                                                    if (seriesNum >= 0) pImage->header_.image_series_index = seriesNum;
                                                }

                                                // set the image data
                                                size_t RO = pImage->get_size(0);
                                                size_t E1 = pImage->get_size(1);
                                                size_t E2 = pImage->get_size(2);

                                                dim3D[0] = RO;
                                                dim3D[1] = E1;
                                                dim3D[2] = E2;

                                                pImage->header_.matrix_size[0] = RO;
                                                pImage->header_.matrix_size[1] = E1;
                                                pImage->header_.matrix_size[2] = E2;

                                                this->decorateImageHeader(pImage->header_, pImage->attrib_, seriesNum, processStr, dataRole, windowCenter, windowWidth, resetImageCommentsParametricMaps, slc, SLC);
                                            }
                                            catch (...)
                                            {
                                                throw;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            Gadgetron::GadgetContainerMessage<hoNDObjectArray< hoMRImage<ValueType, 2> > >* cm1 = new Gadgetron::GadgetContainerMessage<hoNDObjectArray< hoMRImage<ValueType, 2> > >();
            cm1->getObjectPtr()->copyFrom(images);

            if (anchor != NULL)
            {
                if (anchor->putq(cm1) < 0)
                {
                    cm1->release();
                    return false;
                }
            }
            else
            {
                if (this->next()->putq(cm1) < 0)
                {
                    cm1->release();
                    return false;
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadget::sendOutImageBuffer(images, seriesNum, processStr, dataRole) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadget::sendOutImages(Image3DBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole,
        const std::vector<float>& windowCenter, const std::vector<float>& windowWidth, bool resetImageCommentsParametricMaps, Gadget* anchor)
    {
        try
        {
            size_t CHA = images.get_size(0);
            size_t SLC = images.get_size(1);
            size_t CON = images.get_size(2);
            size_t PHS = images.get_size(3);
            size_t REP = images.get_size(4);
            size_t SET = images.get_size(5);
            size_t AVE = images.get_size(6);

            std::string dataRoleString;
            if (!dataRole.empty())
            {
                std::ostringstream ostr;
                for (size_t n = 0; n < dataRole.size(); n++)
                {
                    ostr << dataRole[n] << " - ";
                }

                dataRoleString = ostr.str();
            }

            GDEBUG_CONDITION_STREAM(verbose.value(), "--> GenericImageReconGadget, sending out " << dataRoleString << " images for series " << seriesNum << ", array boundary [CHA SLC CON PHS REP SET AVE] = ["
                << CHA << " " << SLC << " " << CON << " " << PHS << " " << REP << " " << SET << " " << AVE << "] ");

            size_t ave(0), set(0), rep(0), phs(0), con(0), slc(0), cha(0);
            std::vector<size_t> dim3D(3);

            for (ave = 0; ave < AVE; ave++)
            {
                for (set = 0; set < SET; set++)
                {
                    for (rep = 0; rep < REP; rep++)
                    {
                        for (phs = 0; phs < PHS; phs++)
                        {
                            for (con = 0; con < CON; con++)
                            {
                                for (slc = 0; slc < SLC; slc++)
                                {
                                    for (cha = 0; cha < CHA; cha++)
                                    {
                                        hoMRImage<ValueType, 3>* pImage = &images(cha, slc, con, phs, rep, set, ave);
                                        if (pImage != NULL && pImage->get_number_of_elements()>0)
                                        {
                                            T v = Gadgetron::nrm2(*pImage);
                                            if (v < FLT_EPSILON) continue; // do not send out empty image

                                            Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>();
                                            Gadgetron::GadgetContainerMessage<ImgArrayType>* cm2 = new Gadgetron::GadgetContainerMessage<ImgArrayType>();
                                            Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>* cm3 = new Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>();

                                            try
                                            {
                                                cm1->cont(cm2);
                                                cm2->cont(cm3);

                                                // set the ISMRMRD image header
                                                memcpy(cm1->getObjectPtr(), &pImage->header_, sizeof(ISMRMRD::ISMRMRD_ImageHeader));

                                                long long imageNum(0);
                                                if (pImage->attrib_.length(GADGETRON_IMAGENUMBER) == 0)
                                                {
                                                    imageNum = this->computeSeriesImageNumber(*cm1->getObjectPtr(), CHA, cha, 1, 0);
                                                    cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                    pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                }
                                                else
                                                {
                                                    imageNum = pImage->attrib_.as_long(GADGETRON_IMAGENUMBER);
                                                    if (imageNum > 0)
                                                    {
                                                        cm1->getObjectPtr()->image_index = imageNum;
                                                    }
                                                    else
                                                    {
                                                        imageNum = this->computeSeriesImageNumber(*cm1->getObjectPtr(), CHA, cha, 1, 0);
                                                        cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                        pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                    }
                                                }

                                                cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                                                if (add_original_series_num.value())
                                                {
                                                    if (seriesNum >= 0) cm1->getObjectPtr()->image_series_index += seriesNum;
                                                }
                                                else
                                                {
                                                    if (seriesNum >= 0) cm1->getObjectPtr()->image_series_index = seriesNum;
                                                }

                                                // set the image data
                                                size_t RO = pImage->get_size(0);
                                                size_t E1 = pImage->get_size(1);
                                                size_t E2 = pImage->get_size(2);

                                                dim3D[0] = RO;
                                                dim3D[1] = E1;
                                                dim3D[2] = E2;

                                                cm1->getObjectPtr()->matrix_size[0] = RO;
                                                cm1->getObjectPtr()->matrix_size[1] = E1;
                                                cm1->getObjectPtr()->matrix_size[2] = E2;

                                                cm2->getObjectPtr()->create(dim3D);
                                                memcpy(cm2->getObjectPtr()->get_data_ptr(), pImage->get_data_ptr(), pImage->get_number_of_bytes());

                                                // set the attributes
                                                *cm3->getObjectPtr() = pImage->attrib_;

                                                this->decorateImageHeader(*cm1->getObjectPtr(), *cm3->getObjectPtr(), seriesNum, processStr, dataRole, windowCenter, windowWidth, resetImageCommentsParametricMaps, slc, SLC);

                                                if (anchor != NULL)
                                                {
                                                    if (anchor->putq(cm1) < 0)
                                                    {
                                                        cm1->release();
                                                        return false;
                                                    }
                                                }
                                                else
                                                {
                                                    if (this->next()->putq(cm1) < 0)
                                                    {
                                                        cm1->release();
                                                        return false;
                                                    }
                                                }
                                            }
                                            catch (...)
                                            {
                                                cm1->release();
                                                throw;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadget::sendOutImages(images, seriesNum, processStr, dataRole) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadget::sendOutImageBuffer(Image3DBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole,
        const std::vector<float>& windowCenter, const std::vector<float>& windowWidth, bool resetImageCommentsParametricMaps, Gadget* anchor)
    {
        try
        {
            size_t CHA = images.get_size(0);
            size_t SLC = images.get_size(1);
            size_t CON = images.get_size(2);
            size_t PHS = images.get_size(3);
            size_t REP = images.get_size(4);
            size_t SET = images.get_size(5);
            size_t AVE = images.get_size(6);

            std::string dataRoleString;
            if (!dataRole.empty())
            {
                std::ostringstream ostr;
                for (size_t n = 0; n < dataRole.size(); n++)
                {
                    ostr << dataRole[n] << " - ";
                }

                dataRoleString = ostr.str();
            }

            GDEBUG_CONDITION_STREAM(verbose.value(), "--> GenericImageReconGadget, sending out image array " << dataRoleString << " images for series " << seriesNum << ", array boundary [CHA SLC CON PHS REP SET AVE] = ["
                << CHA << " " << SLC << " " << CON << " " << PHS << " " << REP << " " << SET << " " << AVE << "] ");

            size_t ave(0), set(0), rep(0), phs(0), con(0), slc(0), cha(0);
            std::vector<size_t> dim3D(3);

            for (ave = 0; ave < AVE; ave++)
            {
                for (set = 0; set < SET; set++)
                {
                    for (rep = 0; rep < REP; rep++)
                    {
                        for (phs = 0; phs < PHS; phs++)
                        {
                            for (con = 0; con < CON; con++)
                            {
                                for (slc = 0; slc < SLC; slc++)
                                {
                                    for (cha = 0; cha < CHA; cha++)
                                    {
                                        hoMRImage<ValueType, 3>* pImage = &images(cha, slc, con, phs, rep, set, ave);
                                        if (pImage != NULL && pImage->get_number_of_elements() > 0)
                                        {
                                            T v = Gadgetron::nrm2(*pImage);
                                            if (v < FLT_EPSILON) continue; // do not send out empty image

                                            try
                                            {
                                                long long imageNum(0);
                                                if (pImage->attrib_.length(GADGETRON_IMAGENUMBER) == 0)
                                                {
                                                    imageNum = this->computeSeriesImageNumber(pImage->header_, CHA, cha, 1, 0);
                                                    pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                }
                                                else
                                                {
                                                    imageNum = pImage->attrib_.as_long(GADGETRON_IMAGENUMBER);
                                                    if (imageNum <= 0)
                                                    {
                                                        imageNum = this->computeSeriesImageNumber(pImage->header_, CHA, cha, 1, 0);
                                                        pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                    }
                                                }

                                                pImage->header_.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                                                if (add_original_series_num.value())
                                                {
                                                    if (seriesNum >= 0) pImage->header_.image_series_index += seriesNum;
                                                }
                                                else
                                                {
                                                    if (seriesNum >= 0) pImage->header_.image_series_index = seriesNum;
                                                }

                                                // set the image data
                                                size_t RO = pImage->get_size(0);
                                                size_t E1 = pImage->get_size(1);
                                                size_t E2 = pImage->get_size(2);

                                                dim3D[0] = RO;
                                                dim3D[1] = E1;
                                                dim3D[2] = E2;

                                                pImage->header_.matrix_size[0] = RO;
                                                pImage->header_.matrix_size[1] = E1;
                                                pImage->header_.matrix_size[2] = E2;

                                                this->decorateImageHeader(pImage->header_, pImage->attrib_, seriesNum, processStr, dataRole, windowCenter, windowWidth, resetImageCommentsParametricMaps, slc, SLC);
                                            }
                                            catch (...)
                                            {
                                                throw;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            Gadgetron::GadgetContainerMessage<hoNDObjectArray< hoMRImage<ValueType, 3> > >* cm1 = new Gadgetron::GadgetContainerMessage<hoNDObjectArray< hoMRImage<ValueType, 3> > >();
            cm1->getObjectPtr()->copyFrom(images);

            if (anchor != NULL)
            {
                if (anchor->putq(cm1) < 0)
                {
                    cm1->release();
                    return false;
                }
            }
            else
            {
                if (this->next()->putq(cm1) < 0)
                {
                    cm1->release();
                    return false;
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadget::sendOutImageBuffer(images, seriesNum, processStr, dataRole) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadget::sendOutRGBImages(Image3DMagBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, bool resetImageCommentsParametricMaps, Gadget* anchor)
    {
        try
        {
            size_t CHA = images.get_size(0);
            size_t SLC = images.get_size(1);
            size_t CON = images.get_size(2);
            size_t PHS = images.get_size(3);
            size_t REP = images.get_size(4);
            size_t SET = images.get_size(5);
            size_t AVE = images.get_size(6);

            std::string dataRoleString;
            if (!dataRole.empty())
            {
                std::ostringstream ostr;
                for (size_t n = 0; n < dataRole.size(); n++)
                {
                    ostr << dataRole[n] << " - ";
                }

                dataRoleString = ostr.str();
            }

            GDEBUG_CONDITION_STREAM(verbose.value(), "--> GenericImageReconGadget, sending out RBG " << dataRoleString << " images for series " << seriesNum << ", array boundary [CHA SLC CON PHS REP SET AVE] = ["
                << CHA << " " << SLC << " " << CON << " " << PHS << " " << REP << " " << SET << " " << AVE << "] ");

            size_t ave(0), set(0), rep(0), phs(0), con(0), slc(0), cha(0);
            std::vector<size_t> dim3D(3);

            for (ave = 0; ave < AVE; ave++)
            {
                for (set = 0; set < SET; set++)
                {
                    for (rep = 0; rep < REP; rep++)
                    {
                        for (phs = 0; phs < PHS; phs++)
                        {
                            for (con = 0; con < CON; con++)
                            {
                                for (slc = 0; slc < SLC; slc++)
                                {
                                    for (cha = 0; cha < CHA; cha++)
                                    {
                                        hoMRImage<T, 3>* pImage = &images(cha, slc, con, phs, rep, set, ave);
                                        if (pImage != NULL && pImage->get_number_of_elements() > 0)
                                        {
                                            T v = Gadgetron::nrm2(*pImage);
                                            if (v < FLT_EPSILON) continue; // do not send out empty image

                                            Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>();
                                            Gadgetron::GadgetContainerMessage<ImgRGBType>* cm2 = new Gadgetron::GadgetContainerMessage<ImgRGBType>();
                                            Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>* cm3 = new Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>();

                                            try
                                            {
                                                cm1->cont(cm2);
                                                cm2->cont(cm3);

                                                // set the ISMRMRD image header
                                                memcpy(cm1->getObjectPtr(), &pImage->header_, sizeof(ISMRMRD::ISMRMRD_ImageHeader));

                                                long long imageNum(0);
                                                if (pImage->attrib_.length(GADGETRON_IMAGENUMBER) == 0)
                                                {
                                                    imageNum = this->computeSeriesImageNumber(*cm1->getObjectPtr(), CHA, cha, 1, 0);
                                                    cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                    pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                }
                                                else
                                                {
                                                    imageNum = pImage->attrib_.as_long(GADGETRON_IMAGENUMBER);
                                                    if (imageNum > 0)
                                                    {
                                                        cm1->getObjectPtr()->image_index = imageNum;
                                                    }
                                                    else
                                                    {
                                                        imageNum = this->computeSeriesImageNumber(*cm1->getObjectPtr(), CHA, cha, 1, 0);
                                                        cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                        pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                    }
                                                }

                                                cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_USHORT;
                                                if (add_original_series_num.value())
                                                {
                                                    if (seriesNum >= 0) cm1->getObjectPtr()->image_series_index += seriesNum;
                                                }
                                                else
                                                {
                                                    if (seriesNum >= 0) cm1->getObjectPtr()->image_series_index = seriesNum;
                                                }

                                                // set the image data
                                                size_t RO = pImage->get_size(0);
                                                size_t E1 = pImage->get_size(1);
                                                size_t E2 = pImage->get_size(2);

                                                dim3D[0] = RO;
                                                dim3D[1] = E1;
                                                dim3D[2] = E2;

                                                cm1->getObjectPtr()->matrix_size[0] = RO;
                                                cm1->getObjectPtr()->matrix_size[1] = E1;
                                                cm1->getObjectPtr()->matrix_size[2] = E2;

                                                cm2->getObjectPtr()->create(dim3D);
                                                uint16_t* p_img = cm2->getObjectPtr()->get_data_ptr();

                                                for (size_t k = 0; k < RO * E1 * E2; k++)
                                                {
                                                    p_img[k] = uint16_t((*pImage)(k));
                                                }

                                                // set the attributes
                                                *cm3->getObjectPtr() = pImage->attrib_;

                                                std::vector<float> windowCenter, windowWidth;
                                                this->decorateImageHeader(*cm1->getObjectPtr(), *cm3->getObjectPtr(), seriesNum, processStr, dataRole, windowCenter, windowWidth, resetImageCommentsParametricMaps, slc, SLC);

                                                if (anchor != NULL)
                                                {
                                                    if (anchor->putq(cm1) < 0)
                                                    {
                                                        cm1->release();
                                                        return false;
                                                    }
                                                }
                                                else
                                                {
                                                    if (this->next()->putq(cm1) < 0)
                                                    {
                                                        cm1->release();
                                                        return false;
                                                    }
                                                }
                                            }
                                            catch (...)
                                            {
                                                cm1->release();
                                                throw;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadget::sendOutImages(images, seriesNum, processStr, dataRole) ... ");
            return false;
        }

        return true;
    }


    bool GenericImageReconGadget::sendOutRGBImageBuffer(Image3DMagBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, bool resetImageCommentsParametricMaps, Gadget* anchor)
    {
        try
        {
            size_t CHA = images.get_size(0);
            size_t SLC = images.get_size(1);
            size_t CON = images.get_size(2);
            size_t PHS = images.get_size(3);
            size_t REP = images.get_size(4);
            size_t SET = images.get_size(5);
            size_t AVE = images.get_size(6);

            std::string dataRoleString;
            if (!dataRole.empty())
            {
                std::ostringstream ostr;
                for (size_t n = 0; n < dataRole.size(); n++)
                {
                    ostr << dataRole[n] << " - ";
                }

                dataRoleString = ostr.str();
            }

            GDEBUG_CONDITION_STREAM(verbose.value(), "--> GenericImageReconGadget, sending out RGB image array " << dataRoleString << " images for series " << seriesNum << ", array boundary [CHA SLC CON PHS REP SET AVE] = ["
                << CHA << " " << SLC << " " << CON << " " << PHS << " " << REP << " " << SET << " " << AVE << "] ");

            size_t ave(0), set(0), rep(0), phs(0), con(0), slc(0), cha(0);
            std::vector<size_t> dim3D(3);

            for (ave = 0; ave < AVE; ave++)
            {
                for (set = 0; set < SET; set++)
                {
                    for (rep = 0; rep < REP; rep++)
                    {
                        for (phs = 0; phs < PHS; phs++)
                        {
                            for (con = 0; con < CON; con++)
                            {
                                for (slc = 0; slc < SLC; slc++)
                                {
                                    for (cha = 0; cha < CHA; cha++)
                                    {
                                        hoMRImage<T, 3>* pImage = &images(cha, slc, con, phs, rep, set, ave);
                                        if (pImage != NULL && pImage->get_number_of_elements() > 0)
                                        {
                                            T v = Gadgetron::nrm2(*pImage);
                                            if (v < FLT_EPSILON) continue; // do not send out empty image

                                            try
                                            {
                                                long long imageNum(0);
                                                if (pImage->attrib_.length(GADGETRON_IMAGENUMBER) == 0)
                                                {
                                                    imageNum = this->computeSeriesImageNumber(pImage->header_, CHA, cha, 1, 0);
                                                    pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                }
                                                else
                                                {
                                                    imageNum = pImage->attrib_.as_long(GADGETRON_IMAGENUMBER);
                                                    if (imageNum <= 0)
                                                    {
                                                        imageNum = this->computeSeriesImageNumber(pImage->header_, CHA, cha, 1, 0);
                                                        pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                    }
                                                }

                                                pImage->header_.data_type = ISMRMRD::ISMRMRD_USHORT;
                                                if (add_original_series_num.value())
                                                {
                                                    if (seriesNum >= 0) pImage->header_.image_series_index += seriesNum;
                                                }
                                                else
                                                {
                                                    if (seriesNum >= 0) pImage->header_.image_series_index = seriesNum;
                                                }

                                                // set the image data
                                                size_t RO = pImage->get_size(0);
                                                size_t E1 = pImage->get_size(1);
                                                size_t E2 = pImage->get_size(2);

                                                dim3D[0] = RO;
                                                dim3D[1] = E1;
                                                dim3D[2] = E2;

                                                pImage->header_.matrix_size[0] = RO;
                                                pImage->header_.matrix_size[1] = E1;
                                                pImage->header_.matrix_size[2] = E2;

                                                std::vector<float> wc, ww;
                                                this->decorateImageHeader(pImage->header_, pImage->attrib_, seriesNum, processStr, dataRole, wc, ww, resetImageCommentsParametricMaps, slc, SLC);
                                            }
                                            catch (...)
                                            {
                                                throw;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            Gadgetron::GadgetContainerMessage<hoNDObjectArray< hoMRImage<T, 3> > >* cm1 = new Gadgetron::GadgetContainerMessage<hoNDObjectArray< hoMRImage<T, 3> > >();
            cm1->getObjectPtr()->copyFrom(images);

            if (anchor != NULL)
            {
                if (anchor->putq(cm1) < 0)
                {
                    cm1->release();
                    return false;
                }
            }
            else
            {
                if (this->next()->putq(cm1) < 0)
                {
                    cm1->release();
                    return false;
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadget::sendOutRGBImageBuffer(images, seriesNum, processStr, dataRole) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadget::releaseImageBuffer(hoNDObjectArray< hoMRImage<ValueType, 2> >& buf)
    {
        try
        {
            buf.clear();
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadget::releaseImageBuffer(buf) ... ");
            return false;
        }

        return true;
    }

    bool GenericImageReconGadget::releaseImageBuffer(hoNDObjectArray< hoMRImage<ValueType, 3> >& buf)
    {
        try
        {
            buf.clear();
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in GenericImageReconGadget::releaseImageBuffer(buf) ... ");
            return false;
        }

        return true;
    }

    GADGET_FACTORY_DECLARE(GenericImageReconGadget)
}
