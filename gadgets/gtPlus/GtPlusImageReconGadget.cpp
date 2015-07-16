
#include "GtPlusImageReconGadget.h"
#include "GtPlusGadgetOpenMP.h"
#include <iomanip>

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

    GtPlusImageReconGadget::GtPlusImageReconGadget()
    {
        image_series_num_ = 100;

        debugFolder_ = "DebugOutput";

        performTiming_ = true;

        verboseMode_ = false;

        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);

        Gadgetron::prepOpenMP();
    }

    GtPlusImageReconGadget::~GtPlusImageReconGadget()
    {

    }

    bool GtPlusImageReconGadget::readParameters()
    {
        try
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "------> GtPlusImageReconGadget parameters <------");

            verboseMode_ = verboseMode.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "verboseMode_ is " << verboseMode_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            debugFolder_ = debugFolder.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "debugFolder_ is " << debugFolder_);

            if ( !debugFolder_.empty() )
            {
                Gadgetron::getDebugFolderPath(debugFolder_, debugFolder_fullPath_, verboseMode_);
            }
            else
            {
                GDEBUG_STREAM("GtPlusImageRecon, debugFolder is not set ...");
            }

            performTiming_ = performTiming.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "performTiming_ is " << performTiming_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");
        }
        catch(...)
        {
            GERROR_STREAM("Errors in GtPlusImageReconGadget::readParameters() ... ");
            return false;
        }

        return true;
    }

    int GtPlusImageReconGadget::process_config(ACE_Message_Block* mb)
    {
        // read in parameters from the xml
        GADGET_CHECK_RETURN(this->readParameters(), GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try {
          deserialize(mb->rd_ptr(),h);
        } catch (...) {
          GDEBUG("Error parsing ISMRMRD Header");
          throw;
          return GADGET_FAIL;
        }

        // seq object
        if (h.encoding.size() != 1)
        {
            GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
            GDEBUG("This simple GtPlusAccumulatorImageTriggerGadget only supports one encoding space\n");
            return GADGET_FAIL;
        }

        GADGET_CHECK_RETURN(findEncodingLimits(h, meas_max_idx_, verboseMode_), GADGET_FAIL);

        return GADGET_OK;
    }

    int GtPlusImageReconGadget::process(GadgetContainerMessage<ImageBufferType>* m1)
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusImageReconGadget::process(...) starts ... ");

        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        ImageBufferType& ori = *m1->getObjectPtr();

        if ( ori.get_number_of_elements() == 1 )
        {
            size_t num = (*ori(0)).attrib_.length(GADGETRON_DATA_ROLE);
            GADGET_CHECK_RETURN(num>0, GADGET_FAIL);

            dataRole.resize(num);

            for ( size_t ii=0; ii<num; ii++ )
            {
                dataRole[ii] = std::string( (*ori(0)).attrib_.as_str(GADGETRON_DATA_ROLE, ii) );
            }

            if ( (dataRole[0] == GADGETRON_IMAGE_GFACTOR)
                || (dataRole[0] == GADGETRON_IMAGE_SNR_MAP)
                || (dataRole[0] == GADGETRON_IMAGE_STD_MAP)
                || (dataRole[0] == GADGETRON_IMAGE_WRAPAROUNDMAP) )
            {
                GADGET_CHECK_RETURN(this->sendOutImages(ori, image_series_num_++, processStr, dataRole), GADGET_FAIL);
                GADGET_CHECK_RETURN(this->releaseImageBuffer(ori), GADGET_FAIL);
                return GADGET_OK;
            }
        }

        this->processImageBuffer(ori);

        this->releaseImageBuffer(ori);

        m1->release();

        return GADGET_OK;
    }

    int GtPlusImageReconGadget::processImageBuffer(ImageBufferType& ori)
    {
        std::vector<std::string> processStr;
        std::vector<std::string> dataRole;

        boost::shared_ptr< std::vector<size_t> > dims = ori.get_dimensions();
        GDEBUG_CONDITION_STREAM(verboseMode_, "[Cha Slice E2 Con Phase Rep Set Ave] = [" << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " "
            << (*dims)[3] << " " << (*dims)[4]  << " " << (*dims)[5] << " "
            << (*dims)[6] << " " << (*dims)[7] << "]");

        this->sendOutImages(ori, image_series_num_++, processStr, dataRole);

        return GADGET_OK;
    }

    bool GtPlusImageReconGadget::fillWithNULL(ImageBufferType& buf)
    {
        try
        {
            size_t N = buf.get_number_of_elements();
            size_t ii;
            for ( ii=0; ii<N; ii++ )
            {
                buf(ii) = NULL;
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in GtPlusImageReconGadget::fillWithNULL(ImageBufferType& buf) ... ");
            return false;
        }

        return true;
    }

    bool GtPlusImageReconGadget::releaseImageBuffer(ImageBufferType& buf)
    {
        try
        {
            size_t N = buf.get_number_of_elements();
            size_t ii;
            for ( ii=0; ii<N; ii++ )
            {
                ImageType* pImage = buf(ii);
                if ( pImage != NULL )
                {
                    delete pImage;
                    buf(ii) = NULL;
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in GtPlusImageReconGadget::releaseImageBuffer(ImageBufferType& buf) ... ");
            return false;
        }

        return true;
    }

    hoNDImage<std::complex<float>, 3>* GtPlusImageReconGadget::getImage3DFromImage2D(ImageBufferType& ori, size_t cha, size_t slc, size_t con, size_t phs, size_t rep, size_t set, size_t ave)
    {
        ImageType* pImage2D = ori(cha, slc, 0, con, phs, rep, set, ave);
        GADGET_CHECK_THROW(pImage2D!=NULL);

        size_t RO = pImage2D->get_size(0);
        size_t E1 = pImage2D->get_size(1);
        size_t E2 = ori.get_size(2);

        Image3DType* pImage3D = new Image3DType(RO, E1, E2);
        GADGET_CHECK_THROW(pImage3D!=NULL);

        pImage3D->attrib_ = pImage2D->attrib_;

        size_t e2;
        for ( e2=0; e2<E2; e2++ )
        {
            pImage2D = ori(cha, slc, e2, con, phs, rep, set, ave);
            GADGET_CHECK_THROW(pImage2D!=NULL);

            memcpy(pImage3D->begin()+e2*RO*E1, pImage2D->begin(), sizeof(ValueType)*RO*E1 );
        }

        return pImage3D;
    }

    bool GtPlusImageReconGadget::getImage2DFromImage3D(Image3DType& image3D, ImageBufferType& image2DBuf)
    {
        size_t RO = image3D.get_size(0);
        size_t E1 = image3D.get_size(1);
        size_t E2 = image3D.get_size(2);

        std::vector<size_t> dim(1);
        dim[0] = E2;
        image2DBuf.create(dim);

        size_t e2;
        for ( e2=0; e2<E2; e2++ )
        {
            ImageType* pImage2D = new ImageType(RO, E1);
            GADGET_CHECK_RETURN_FALSE(pImage2D!=NULL);

            memcpy(pImage2D->begin(), image3D.begin()+e2*RO*E1, sizeof(ValueType)*RO*E1);

            image2DBuf(e2) = pImage2D;
        }

        return true;
    }

    size_t GtPlusImageReconGadget::computeSeriesImageNumber (ISMRMRD::ImageHeader& imheader, size_t nCHA, size_t cha, size_t nE2, size_t e2)
    {
        size_t nSET = meas_max_idx_.set+1;
        size_t nREP = meas_max_idx_.repetition+1;
        size_t nPHS = meas_max_idx_.phase+1;
        size_t nSLC = meas_max_idx_.slice+1;
        size_t nCON = meas_max_idx_.contrast+1;
        if ( nE2 == 0 ) nE2 = 1;

        size_t imageNum = imheader.average*nREP*nSET*nPHS*nCON*nSLC*nE2*nCHA
            + imheader.repetition*nSET*nPHS*nCON*nSLC*nE2*nCHA
            + imheader.set*nPHS*nCON*nSLC*nE2*nCHA
            + imheader.phase*nCON*nSLC*nE2*nCHA
            + imheader.contrast*nSLC*nE2*nCHA
            + imheader.slice*nE2*nCHA
            + e2*nCHA
            + cha
            + 1;

        return imageNum;
    }

    bool GtPlusImageReconGadget::sendOutImages(ImageBufferType& images, int seriesNum, const std::vector<std::string>& processStr, const std::vector<std::string>& dataRole, const std::vector<float>& windowCenter, const std::vector<float>& windowWidth)
    {
        try
        {
            size_t CHA = images.get_size(0);
            size_t SLC = images.get_size(1);
            size_t E2  = images.get_size(2);
            size_t CON = images.get_size(3);
            size_t PHS = images.get_size(4);
            size_t REP = images.get_size(5);
            size_t SET = images.get_size(6);
            size_t AVE = images.get_size(7);

            GDEBUG_CONDITION_STREAM(verboseMode_, "--> GtPlusImageReconGadget, sending out images, array boundary [CHA SLC E2 CON PHS REP SET AVE] = ["
                << CHA << " " << SLC << " "
                << E2 << " " << CON << " "
                << PHS << " " << REP << " "
                << SET << " " << AVE << "] " );

            size_t ave(0), set(0), rep(0), phs(0), con(0), e2(0), slc(0), cha(0);
            std::vector<size_t> dim2D(2);

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
                                for ( e2=0; e2<E2; e2++ )
                                {
                                    for ( slc=0; slc<SLC; slc++ )
                                    {
                                        for ( cha=0; cha<CHA; cha++ )
                                        {
                                            ImageType* pImage = images(cha, slc, e2, con, phs, rep, set, ave);
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
                                                    GADGET_CHECK_THROW( gtPlus_util_.setImageHeaderISMRMRDFromMetaAttributes(pImage->attrib_, *cm1->getObjectPtr()) );

                                                    //long long imageNum(0);
                                                    //if ( pImage->attrib_.attributeInteger_.get(GADGETRON_IMAGENUMBER, 0, imageNum) )
                                                    //{
                                                    //    cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                    //}

                                                    long long imageNum = this->computeSeriesImageNumber (*cm1->getObjectPtr(), CHA, cha, E2, e2);
                                                    cm1->getObjectPtr()->image_index = (uint16_t)imageNum;
                                                    pImage->attrib_.set(GADGETRON_IMAGENUMBER, (long)imageNum);

                                                    cm1->getObjectPtr()->image_series_index = seriesNum;

                                                    // set the image data
                                                    size_t RO = pImage->get_size(0);
                                                    size_t E1 = pImage->get_size(1);

                                                    dim2D[0] = RO;
                                                    dim2D[1] = E1;

                                                    cm2->getObjectPtr()->create(dim2D);
                                                    memcpy(cm2->getObjectPtr()->get_data_ptr(), pImage->get_data_ptr(), pImage->get_number_of_bytes());

                                                    // set the attributes
                                                    *cm3->getObjectPtr() = pImage->attrib_;

                                                    if ( !dataRole.empty() && (dataRole[0]!=GADGETRON_IMAGE_REGULAR) )
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
                                                        bool isParametricT2StarMaskMap = false;
                                                        bool isParametricT2StarSDMap = false;
                                                        bool isParametricT2StarAMap = false;
                                                        bool isParametricT2StarTruncMap = false;

                                                        if ( !dataRole.empty() )
                                                        {
                                                            size_t n;
                                                            for ( n=0; n<dataRole.size(); n++ )
                                                            {
                                                                if ( dataRole[n] == GADGETRON_IMAGE_PSIR )
                                                                {
                                                                    isRealImage = true;
                                                                }

                                                                if ( (dataRole[n]==GADGETRON_IMAGE_T1MAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_T1SDMAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_T2MAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_T2SDMAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_T2STARMAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_T2STARMASKMAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_T2STARSDMAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_T2STARAMAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_T2STARTRUNCMAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_FREQMAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_B1MAP)
                                                                    || (dataRole[n]==GADGETRON_IMAGE_FLIPANGLEMAP) )
                                                                {
                                                                    isParametricMap = true;
                                                                }

                                                                if ( dataRole[n]==GADGETRON_IMAGE_T1MAP )
                                                                {
                                                                    isParametricT1Map = true;
                                                                }

                                                                if ( dataRole[n]==GADGETRON_IMAGE_T1SDMAP )
                                                                {
                                                                    isParametricT1SDMap = true;
                                                                }

                                                                if ( dataRole[n]==GADGETRON_IMAGE_T2MAP )
                                                                {
                                                                    isParametricT2Map = true;
                                                                }

                                                                if ( dataRole[n]==GADGETRON_IMAGE_T2SDMAP )
                                                                {
                                                                    isParametricT2SDMap = true;
                                                                }

                                                                if ( dataRole[n]==GADGETRON_IMAGE_T2STARMAP )
                                                                {
                                                                    isParametricT2StarMap = true;
                                                                }

                                                                if ( dataRole[n]==GADGETRON_IMAGE_T2STARSDMAP )
                                                                {
                                                                    isParametricT2StarSDMap = true;
                                                                }

                                                                if ( dataRole[n]==GADGETRON_IMAGE_T2STARAMAP )
                                                                {
                                                                    isParametricT2StarAMap = true;
                                                                }

                                                                if ( dataRole[n]==GADGETRON_IMAGE_T2STARTRUNCMAP )
                                                                {
                                                                    isParametricT2StarTruncMap = true;
                                                                }

                                                                if ( dataRole[n]==GADGETRON_IMAGE_T2STARMASKMAP )
                                                                {
                                                                    isParametricT2StarMaskMap = true;
                                                                }
                                                            }

                                                            std::vector<std::string> dataRoleAll;
                                                            Gadgetron::getISMRMRMetaValues(*cm3->getObjectPtr(), GADGETRON_DATA_ROLE, dataRoleAll);

                                                            if ( !debugFolder_fullPath_.empty() )
                                                            {
                                                                std::ostringstream ostr;
                                                                for ( n=0; n<dataRoleAll.size(); n++ )
                                                                {
                                                                    ostr << dataRoleAll[n] << "_";
                                                                }
                                                                ostr << cm1->getObjectPtr()->image_index;

                                                                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(*cm2->getObjectPtr(), debugFolder_fullPath_+ostr.str()); }
                                                            }

                                                            // double check the image type
                                                            if ( isRealImage )
                                                            {
                                                                cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_REAL;
                                                            }

                                                            // image comment
                                                            if ( isParametricMap )
                                                            {
                                                                // reset the image comment for maps

                                                                std::vector<std::string> commentStr(dataRole.size()+1);

                                                                commentStr[0] = "GT";
                                                                for ( n=0; n<dataRole.size(); n++ )
                                                                {
                                                                    commentStr[n+1] = dataRole[n];
                                                                }

                                                                Gadgetron::setISMRMRMetaValues(*cm3->getObjectPtr(), GADGETRON_IMAGECOMMENT, commentStr);

                                                                // get the scaling ratio
                                                                float scalingRatio = 1;
                                                                try
                                                                {
                                                                    scalingRatio = (float)(cm3->getObjectPtr()->as_double(GADGETRON_IMAGE_SCALE_RATIO, 0));

                                                                    std::ostringstream ostr;
                                                                    ostr << "x" << scalingRatio;
                                                                    std::string scalingStr = ostr.str();
                                                                    cm3->getObjectPtr()->append(GADGETRON_IMAGECOMMENT, scalingStr.c_str());

                                                                    if ( isParametricT1Map || isParametricT1SDMap || isParametricT2Map || isParametricT2SDMap || isParametricT2StarMap || isParametricT2StarSDMap )
                                                                    {
                                                                        std::ostringstream ostr;
                                                                        ostr << std::setprecision(3) << 1.0f/scalingRatio << "ms";
                                                                        std::string unitStr = ostr.str();

                                                                        cm3->getObjectPtr()->append(GADGETRON_IMAGECOMMENT, unitStr.c_str());
                                                                    }
                                                                }
                                                                catch(...)
                                                                {
                                                                    GWARN_STREAM("Image attrib does not have the scale ratio ...");
                                                                    scalingRatio = 1;
                                                                }
                                                            }
                                                            else
                                                            {
                                                                for ( n=0; n<dataRole.size(); n++ )
                                                                {
                                                                    cm3->getObjectPtr()->append(GADGETRON_IMAGECOMMENT, dataRole[n].c_str());
                                                                }
                                                            }

                                                            // seq description
                                                            Gadgetron::appendISMRMRMetaValues(*cm3->getObjectPtr(), GADGETRON_SEQUENCEDESCRIPTION, dataRoleAll);
                                                        }

                                                        GDEBUG_CONDITION_STREAM(verboseMode_, "--> GtPlusImageReconGadget, sending out 2D image [CHA SLC E2 CON PHS REP SET AVE] = ["
                                                            << cha << " "
                                                            << cm1->getObjectPtr()->slice << " "
                                                            << e2 << " "
                                                            << cm1->getObjectPtr()->contrast << " "
                                                            << cm1->getObjectPtr()->phase << " "
                                                            << cm1->getObjectPtr()->repetition << " "
                                                            << cm1->getObjectPtr()->set << " "
                                                            << cm1->getObjectPtr()->average << "] \t"
                                                            << " -- Image number -- " << cm1->getObjectPtr()->image_index);

                                                        // image processing history
                                                        if ( !processStr.empty() )
                                                        {
                                                            size_t n;
                                                            for ( n=0; n<processStr.size(); n++ )
                                                            {
                                                                cm3->getObjectPtr()->append(GADGETRON_IMAGEPROCESSINGHISTORY, processStr[n].c_str());
                                                            }
                                                        }

                                                        if ( windowCenter.size()==SLC && windowWidth.size()==SLC )
                                                        {
                                                            cm3->getObjectPtr()->set(GADGETRON_IMAGE_WINDOWCENTER, (long)windowCenter[slc]);
                                                            cm3->getObjectPtr()->set(GADGETRON_IMAGE_WINDOWWIDTH, (long)windowWidth[slc]);
                                                        }
                                                    }

                                                    if ( this->next()->putq(cm1) < 0 )
                                                    {
                                                        cm1->release();
                                                        return false;
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
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in GtPlusImageReconGadget::sendOutImages(images, seriesNum, processStr, dataRole) ... ");
            return false;
        }

        return true;
    }

    int GtPlusImageReconGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GtPlusImageReconGadget - close(flags) : " << flags);

        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

        if ( flags != 0 )
        {
            std::string procTime;
            gtPlus_util_.getCurrentMoment(procTime);

            GDEBUG_STREAM("* ============================================================================== *");
            GDEBUG_STREAM("---> Image recon phase, Current processing time : " << procTime << " <---");
            GDEBUG_STREAM("* ============================================================================== *");
        }

        return GADGET_OK;
    }

    bool GtPlusImageReconGadget::exportImageContainer2D(ImageContainer2DType& input, const std::string& prefix)
    {
        if ( !this->debugFolder_.empty() )
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<ValueType> outArray;

            for ( r=0; r<R; r++ )
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if ( !this->debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(outArray, this->debugFolder_fullPath_+ostr.str()); }
            }
        }

        return true;
    }

    bool GtPlusImageReconGadget::exportImageContainer2D(ImageContainer2DMagType& input, const std::string& prefix)
    {
        if ( !this->debugFolder_.empty() )
        {
            size_t R = input.rows();

            size_t r;

            hoNDArray<T> outArray;

            for ( r=0; r<R; r++ )
            {
                input.to_NDArray(r, outArray);

                std::ostringstream ostr;
                ostr << prefix << "_" << r;

                if ( !this->debugFolder_fullPath_.empty() ) { gt_exporter_.exportArray(outArray, this->debugFolder_fullPath_+ostr.str()); }
            }
        }

        return true;
    }

    bool GtPlusImageReconGadget::exportImageContainer3D(ImageContainer3DType& input, const std::string& prefix)
    {
        if ( !this->debugFolder_.empty() )
        {
            size_t R = input.rows();

            size_t r, c;
            for ( r=0; r<R; r++ )
            {
                for ( c=0; c<input.cols(r); c++ )
                {
                    std::ostringstream ostr;
                    ostr << prefix << "_" << r << "_" << c;

                    if ( !this->debugFolder_fullPath_.empty() )
                    {
                        gt_exporter_.exportImageComplex(input(r, c), this->debugFolder_fullPath_+ostr.str());
                    }
                }
            }
        }

        return true;
    }

    bool GtPlusImageReconGadget::exportImageContainer3D(ImageContainer3DMagType& input, const std::string& prefix)
    {
        if ( !this->debugFolder_.empty() )
        {
            size_t R = input.rows();

            size_t r, c;
            for ( r=0; r<R; r++ )
            {
                for ( c=0; c<input.cols(r); c++ )
                {
                    std::ostringstream ostr;
                    ostr << prefix << "_" << r << "_" << c;

                    gt_exporter_.exportImage(input(r, c), this->debugFolder_fullPath_+ostr.str());
                }
            }
        }

        return true;
    }

}
