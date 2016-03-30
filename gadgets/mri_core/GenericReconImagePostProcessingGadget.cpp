
#include "GenericReconImagePostProcessingGadget.h"

namespace Gadgetron { 

    template<typename T, int D> 
    GenericReconImagePostProcessingGadget<T, D>::GenericReconImagePostProcessingGadget() : BaseClass(), process_called_times_(0)
    {
        gt_timer_.set_timing_in_destruction(false);
        gt_timer_local_.set_timing_in_destruction(false);
    }

    template<typename T, int D>
    GenericReconImagePostProcessingGadget<T, D>::~GenericReconImagePostProcessingGadget()
    {

    }

    template<typename T, int D> 
    int GenericReconImagePostProcessingGadget<T, D>::process_config(ACE_Message_Block* mb)
    {
        ISMRMRD::IsmrmrdHeader h;
        try
        {
          deserialize(mb->rd_ptr(),h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
            throw;
            return GADGET_FAIL;
        }

        ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
        ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
        ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

        meas_max_idx_.kspace_encode_step_1 = (uint16_t)e_space.matrixSize.y - 1;

        meas_max_idx_.set = (e_limits.set && (e_limits.set->maximum>0)) ? e_limits.set->maximum : 0;
        meas_max_idx_.phase = (e_limits.phase && (e_limits.phase->maximum>0)) ? e_limits.phase->maximum : 0;

        meas_max_idx_.kspace_encode_step_2 = (uint16_t)e_space.matrixSize.z - 1;

        meas_max_idx_.contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum : 0;

        meas_max_idx_.slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;

        meas_max_idx_.repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;

        meas_max_idx_.average = e_limits.average ? e_limits.average->maximum : 0;

        meas_max_idx_.segment = 0;

        // ------------------------------------------------------------

        /*if (!debug_folder.value().empty())
        {
            Gadgetron::get_debug_folder_path(debug_folder.value(), debug_folder_full_path_);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is " << debug_folder_full_path_);
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is not set ... ");
        }*/

        return GADGET_OK;
    }

    template<typename T, int D> 
    int GenericReconImagePostProcessingGadget<T, D>::process(GadgetContainerMessage<ImageBufferType>* m1)
    {
        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconImagePostProcessingGadget::process(...) starts ... ");

        process_called_times_++;

        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        ImageBufferType& ori = *m1->getObjectPtr();

        if ( ori.get_number_of_elements() > 0 )
        {
            size_t num = (*ori(0)).attrib_.length(GADGETRON_DATA_ROLE);
            GADGET_CHECK_RETURN(num>0, GADGET_FAIL);

            dataRole.resize(num);

            for ( size_t ii=0; ii<num; ii++ )
            {
                dataRole[ii] = std::string( (*ori(0)).attrib_.as_str(GADGETRON_DATA_ROLE, ii) );
            }

            if ( dataRole[0] == GADGETRON_IMAGE_GFACTOR )
            {
                GADGET_CHECK_RETURN( this->store_gfactor_map( (*ori(0)) ), GADGET_FAIL);

                if ( send_out_gfactor_map.value() )
                {
                    GDEBUG_STREAM("Image Recon, send out received gfactor map ...");

                    dataRole.clear();
                    GADGET_CHECK_RETURN(this->send_out_images(ori, ori(0)->header_.image_series_index, processStr, dataRole)==GADGET_OK, GADGET_FAIL);
                    GADGET_CHECK_RETURN(this->release_image_buffer(ori) == GADGET_OK, GADGET_FAIL);
                }

                return GADGET_OK;
            }
            else if ( dataRole[0] == GADGETRON_IMAGE_SNR_MAP )
            {
                if ( send_out_snr_map.value() )
                {
                    dataRole.clear();
                    GADGET_CHECK_RETURN(this->send_out_images(ori, ori(0)->header_.image_series_index, processStr, dataRole) == GADGET_OK, GADGET_FAIL);
                    GADGET_CHECK_RETURN(this->release_image_buffer(ori) == GADGET_OK, GADGET_FAIL);
                }

                return GADGET_OK;
            }
            else if ( dataRole[0] == GADGETRON_IMAGE_STD_MAP )
            {
                if ( send_out_std_map.value())
                {
                    dataRole.clear();
                    GADGET_CHECK_RETURN(this->send_out_images(ori, ori(0)->header_.image_series_index, processStr, dataRole) == GADGET_OK, GADGET_FAIL);
                    GADGET_CHECK_RETURN(this->release_image_buffer(ori) == GADGET_OK, GADGET_FAIL);
                }

                return GADGET_OK;
            }

            num = (*ori(0)).attrib_.length(GADGETRON_PASS_IMMEDIATE);
            if ( num > 0 )
            {
                long pass_image_immediately = (*ori(0)).attrib_.as_long(GADGETRON_PASS_IMMEDIATE, 0);
                if ( pass_image_immediately )
                {
                    dataRole.clear();
                    GADGET_CHECK_RETURN(this->send_out_images(ori, ori(0)->header_.image_series_index, processStr, dataRole) == GADGET_OK, GADGET_FAIL);
                    GADGET_CHECK_RETURN(this->release_image_buffer(ori) == GADGET_OK, GADGET_FAIL);
                    return GADGET_OK;
                }
            }
        }

        this->process_image_buffer(ori);

        this->release_image_buffer(ori);

        m1->release();

        return GADGET_OK;
    }

    template<typename T, int D> 
    int GenericReconImagePostProcessingGadget<T, D>::process_image_buffer(ImageBufferType& ori)
    {
        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        size_t RO = ori(0)->get_size(0);
        size_t E1 = ori(0)->get_size(1);
        size_t E2 = ori(0)->get_size(2);

        boost::shared_ptr< std::vector<size_t> > dims = ori.get_dimensions();
        GDEBUG_CONDITION_STREAM(verbose.value(), "[RO E1 E2 Cha Slice Con Phase Rep Set Ave] = [" 
                                                << RO << " " << E1 << " " << E2 << " " 
                                                << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " 
                                                << (*dims)[3] << " " << (*dims)[4]  << " " << (*dims)[5] << " " 
                                                << (*dims)[6] << " " << (*dims)[7] << "]");

        this->send_out_images(ori, processing_result_image_series_num.value(), processStr, dataRole);

        return GADGET_OK;
    }

    template<typename T, int D> 
    size_t GenericReconImagePostProcessingGadget<T, D>::compute_image_number(ISMRMRD::ImageHeader& imheader, size_t nCHA, size_t cha)
    {
        size_t nSET = meas_max_idx_.set+1;
        size_t nREP = meas_max_idx_.repetition+1;
        size_t nPHS = meas_max_idx_.phase+1;
        size_t nSLC = meas_max_idx_.slice+1;
        size_t nCON = meas_max_idx_.contrast+1;

        size_t imageNum = imheader.average*nREP*nSET*nPHS*nCON*nSLC*nCHA
            + imheader.repetition*nSET*nPHS*nCON*nSLC*nCHA
            + imheader.set*nPHS*nCON*nSLC*nCHA
            + imheader.phase*nCON*nSLC*nCHA
            + imheader.contrast*nSLC*nCHA
            + imheader.slice*nCHA
            + cha
            + 1;

        return imageNum;
    }

    template<typename T, int D> 
    int GenericReconImagePostProcessingGadget<T, D>::send_out_images(ImageBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter, const std::vector<float>& windowWidth)
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

            GDEBUG_CONDITION_STREAM(verbose.value(), "--> GenericReconImagePostProcessingGadget, sending out " << dataRoleString << " images for series " << seriesNum << ", array boundary [CHA SLC CON PHS REP SET AVE] = ["
                << CHA << " " << SLC << " " << CON << " " << PHS << " " << REP << " " << SET << " " << AVE << "] ");

            size_t ave(0), set(0), rep(0), phs(0), con(0), slc(0), cha(0);
            std::vector<size_t> dim3D(3);

            for ( ave=0; ave<AVE; ave++ )
            {
                for ( set=0; set<SET; set++ )
                {
                    for ( rep=0; rep<REP; rep++ )
                    {
                        for ( phs=0; phs<PHS; phs++ )
                        {
                            for ( con=0; con<CON; con++ )
                            {
                                for ( slc=0; slc<SLC; slc++ )
                                {
                                    for ( cha=0; cha<CHA; cha++ )
                                    {
                                        ImageType* pImage = images(cha, slc, con, phs, rep, set, ave);
                                        if ( pImage != NULL )
                                        {
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
                                                    imageNum = this->compute_image_number(*cm1->getObjectPtr(), CHA, cha);
                                                    cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                    pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                }
                                                else
                                                {
                                                    imageNum = pImage->attrib_.as_long(GADGETRON_IMAGENUMBER);
                                                    if (imageNum > 0)
                                                    {
                                                        cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                    }
                                                    else
                                                    {
                                                        imageNum = this->compute_image_number(*cm1->getObjectPtr(), CHA, cha);
                                                        cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                        pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);
                                                    }
                                                }

                                                cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                                                cm1->getObjectPtr()->image_series_index = seriesNum;

                                                // set the image data
                                                size_t RO = pImage->get_size(0);
                                                size_t E1 = pImage->get_size(1);
                                                size_t E2 = pImage->get_size(2);

                                                dim3D[0] = RO;
                                                dim3D[1] = E1;
                                                dim3D[2] = E2;

                                                cm1->getObjectPtr()->matrix_size[0] = (uint16_t)RO;
                                                cm1->getObjectPtr()->matrix_size[1] = (uint16_t)E1;
                                                cm1->getObjectPtr()->matrix_size[2] = (uint16_t)E2;

                                                cm2->getObjectPtr()->create(dim3D);
                                                memcpy(cm2->getObjectPtr()->get_data_ptr(), pImage->get_data_ptr(), pImage->get_number_of_bytes());

                                                // set the attributes
                                                *cm3->getObjectPtr() = pImage->attrib_;

                                                if ( this->next()->putq(cm1) < 0 ) 
                                                {
                                                    cm1->release();
                                                    return GADGET_FAIL;
                                                }
                                            }
                                            catch(...)
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
        catch(...)
        {
            GERROR_STREAM("Errors happened in GenericReconImagePostProcessingGadget::sendOutImages(images, seriesNum, processStr, dataRole) ... ");
            return false;
        }

        return GADGET_OK;
    }

    template<typename T, int D> 
    int GenericReconImagePostProcessingGadget<T, D>::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GenericReconImagePostProcessingGadget - close(flags) : " << flags);

        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

        if ( flags != 0 )
        {
        }

        return GADGET_OK;
    }

    template<typename T, int D> 
    int GenericReconImagePostProcessingGadget<T, D>::store_gfactor_map(const ImageType& gmap)
    {
        try
        {
            if ( gfactor_buf_.empty() )
            {
                gfactor_buf_.push_back(gmap);
                return GADGET_OK;
            }

            size_t N = gfactor_buf_.size();

            size_t n;
            for ( n=0; n<N; n++ )
            {
                long nSLC = gfactor_buf_[n].header_.slice;
                long slc = gmap.header_.slice;

                if ( nSLC==slc )
                {
                    return GADGET_OK;
                }
            }

            gfactor_buf_.push_back(gmap);
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in GenericReconImagePostProcessingGadget::store_gfactor_map(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    template<typename T, int D> 
    int GenericReconImagePostProcessingGadget<T, D>::find_gfactor_map(size_t slc, ImageType& gmap)
    {
        if ( gfactor_buf_.empty() ) return GADGET_FAIL;

        size_t N = gfactor_buf_.size();

        size_t n;
        for ( n=0; n<N; n++ )
        {
            long nSLC = gfactor_buf_[n].header_.slice;

            if ( nSLC==slc )
            {
                gmap = gfactor_buf_[n];
                return GADGET_OK;
            }
        }

        return GADGET_FAIL;
    }

    /*template<typename T, int D>
    void GenericReconImagePostProcessingGadget<T, D>::export_image_container(const hoNDImageContainer2D< hoNDImage<std::complex<T>, D> >& input, const std::string& name)
    {
        if (!this->debugFolder_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<ValueType> outArray;

            for (r = 0; r<R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debugFolder_fullPath_.empty()) gt_exporter_.exportArrayComplex(outArray, debugFolder_fullPath_ + ostr.str());
            }
        }
    }

    template<typename T, int D>
    void GenericReconImagePostProcessingGadget<T, D>::export_image_container(const hoNDImageContainer2D< hoNDImage<T, D> >& input, const std::string& name)
    {
        if (!this->debugFolder_.empty())
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<T> outArray;

            for (r = 0; r<R; r++)
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if (!debugFolder_fullPath_.empty()) gt_exporter_.exportArray(outArray, debugFolder_fullPath_ + ostr.str());
            }
        }
    }*/

    GenericReconImagePostProcessing2DGadget::GenericReconImagePostProcessing2DGadget() : BaseClass()
    {
    }

    GenericReconImagePostProcessing2DGadget::~GenericReconImagePostProcessing2DGadget()
    {
    }

    GenericReconImagePostProcessing3DGadget::GenericReconImagePostProcessing3DGadget() : BaseClass()
    {
    }

    GenericReconImagePostProcessing3DGadget::~GenericReconImagePostProcessing3DGadget()
    {
    }

    GADGET_FACTORY_DECLARE(GenericReconImagePostProcessing2DGadget)
    GADGET_FACTORY_DECLARE(GenericReconImagePostProcessing3DGadget)
}
