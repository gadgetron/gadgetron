
#include "GenericReconAccumulateImageTriggerGadget.h"

namespace Gadgetron { 

    template <typename T, int D> 
    GenericReconAccumulateImageTriggerGadget<T, D>::GenericReconAccumulateImageTriggerGadget() : BaseClass(), image_counter_(0), triggered_in_close_(false)
    {
        cha_trigger_ = false;
        slc_trigger_ = false;
        con_trigger_ = false;
        phs_trigger_ = false;
        rep_trigger_ = false;
        set_trigger_ = false;
        ave_trigger_ = false;

        num_of_dimensions_ = 7; // [CHA SLC CON PHS REP SET AVE]

        pass_image_immediate_ = false;
    }

    template <typename T, int D>
    GenericReconAccumulateImageTriggerGadget<T, D>::~GenericReconAccumulateImageTriggerGadget()
    {

    }

    template <typename T, int D>
    int GenericReconAccumulateImageTriggerGadget<T, D>::process_config(ACE_Message_Block* mb)
    {
        cha_trigger_ = TriggerChannel.value();
        slc_trigger_ = TriggerSlice.value();
        con_trigger_ = TriggerContrast.value();
        phs_trigger_ = TriggerPhase.value();
        rep_trigger_ = TriggerRepetition.value();
        set_trigger_ = TriggerSet.value();
        ave_trigger_ = TriggerAverage.value();

        pass_image_immediate_ = PassImageImmediately.value();

        // ---------------------------------------------------------------------------------------------------------
        // pass the xml file
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

        int num_of_conc = 1;
        int retro_gated_images = -1;
        if (h.userParameters)
        {
            for (std::vector<ISMRMRD::UserParameterLong>::const_iterator i = h.userParameters->userParameterLong.begin(); i != h.userParameters->userParameterLong.end(); ++i)
            {
                if (std::strcmp(i->name.c_str(), "NumOfConcatenations") == 0)
                {
                    num_of_conc = i->value;
                }

                if (std::strcmp(i->name.c_str(), "RetroGatedImages") == 0)
                {
                    retro_gated_images = i->value;
                }
            }
        }

        size_t NE = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        // ---------------------------------------------------------------------------------------------------------

        // encoding limits
        ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
        ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
        ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

        meas_max_idx_.kspace_encode_step_1 = (uint16_t)(h.encoding[0].reconSpace.matrixSize.y);

        meas_max_idx_.set                  = (e_limits.set && (e_limits.set->maximum>0)) ? e_limits.set->maximum : 0;
        meas_max_idx_.phase                = (e_limits.phase && (e_limits.phase->maximum>0)) ? e_limits.phase->maximum : 0;

        meas_max_idx_.kspace_encode_step_2 = (uint16_t)(h.encoding[0].reconSpace.matrixSize.z);

        meas_max_idx_.contrast             = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum : 0;

        meas_max_idx_.slice                = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;

        meas_max_idx_.repetition           = e_limits.repetition ? e_limits.repetition->maximum : 0;

        meas_max_idx_.average              = e_limits.average ? e_limits.average->maximum : 0;

        // always combine the SEG
        meas_max_idx_.segment              = 0;

        // if user suplies the dimension limits, read them in
        int meas_max_kspace_encode_step_1_ = meas_max_kspace_encode_step_1.value();
        if (meas_max_kspace_encode_step_1_ > 0) meas_max_idx_.kspace_encode_step_1 = meas_max_kspace_encode_step_1_;

        int meas_max_set_ = meas_max_set.value();
        if (meas_max_set_ > 0) meas_max_idx_.set = meas_max_set_;

        int meas_max_phase_ = meas_max_phase.value();
        if (meas_max_phase_ > 0) meas_max_idx_.phase = meas_max_phase_;
        if (retro_gated_images>0)
        {
            meas_max_phase_ = retro_gated_images - 1;
            meas_max_idx_.phase = meas_max_phase_;
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Max phase is changed to  " << meas_max_idx_.phase);
        }

        int meas_max_kspace_encode_step_2_ = meas_max_kspace_encode_step_2.value();
        if (meas_max_kspace_encode_step_2_ > 0) meas_max_idx_.kspace_encode_step_2 = meas_max_kspace_encode_step_2_;

        int meas_max_contrast_ = meas_max_contrast.value();
        if (meas_max_contrast_ > 0) meas_max_idx_.contrast = meas_max_contrast_;

        int meas_max_slice_ = meas_max_slice.value();
        if (meas_max_slice_ > 0) meas_max_idx_.slice = meas_max_slice_;

        int meas_max_repetition_ = meas_max_repetition.value();
        if (meas_max_repetition_ > 0) meas_max_idx_.repetition = meas_max_repetition_;

        int meas_max_average_ = meas_max_average.value();
        if (meas_max_average_ > 0) meas_max_idx_.average = meas_max_average_;

        if (concatenation_on_repetition.value())
        {
            meas_max_repetition_ = (meas_max_idx_.repetition + 1)*num_of_conc - 1;
            meas_max_idx_.repetition = meas_max_repetition_;
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Max repetition is changed to  " << meas_max_repetition_);
        }

        // allocate the image buffers
        // [Cha Slice Con Phase Rep Set Ave]
        //   0    1    2   3    4   5   6

        dimensions_.resize(GT_DIM_NUM_IMAGE, 0);
        dimensions_[0] = 1;
        dimensions_[1] = meas_max_idx_.slice+1;
        dimensions_[2] = meas_max_idx_.contrast+1;
        dimensions_[3] = meas_max_idx_.phase+1;
        dimensions_[4] = meas_max_idx_.repetition+1;
        dimensions_[5] = meas_max_idx_.set+1;
        dimensions_[6] = meas_max_idx_.average+1;

        imageBuffer_.clear();
        imageSent_.clear();

        otherBuffer_.clear();
        otherSent_.clear();

        // set the dimensions under/not under trigger
        this->set_dimensions_under_trigger();

        GDEBUG_CONDITION_STREAM(this->verbose.value(), "dimension limits                [Slice Con Phase Rep Set Ave] = [" 
                                   << " " << dimensions_[1] 
                                   << " " << dimensions_[2] 
                                   << " " << dimensions_[3]
                                   << " " << dimensions_[4]
                                   << " " << dimensions_[5]
                                   << " " << dimensions_[6] << "]");

        GDEBUG_CONDITION_STREAM(this->verbose.value(), "dimension under trigger         [Cha Slice Con Phase Rep Set Ave] = [" 
                                   << " " << dim_under_trigger_[0] 
                                   << " " << dim_under_trigger_[1] 
                                   << " " << dim_under_trigger_[2] 
                                   << " " << dim_under_trigger_[3]
                                   << " " << dim_under_trigger_[4]
                                   << " " << dim_under_trigger_[5]
                                   << " " << dim_under_trigger_[6] << "]");

        GDEBUG_CONDITION_STREAM(this->verbose.value(), "dimension limits under trigger  [Slice Con Phase Rep Set Ave] = [" 
                                   << " " << dim_limit_under_trigger_[1] 
                                   << " " << dim_limit_under_trigger_[2] 
                                   << " " << dim_limit_under_trigger_[3]
                                   << " " << dim_limit_under_trigger_[4]
                                   << " " << dim_limit_under_trigger_[5]
                                   << " " << dim_limit_under_trigger_[6] << "]");

        GDEBUG_CONDITION_STREAM(this->verbose.value(), "dimension NOT under trigger     [Cha Slice Con Phase Rep Set Ave] = [" 
                                   << " " << dim_not_under_trigger_[0] 
                                   << " " << dim_not_under_trigger_[1] 
                                   << " " << dim_not_under_trigger_[2] 
                                   << " " << dim_not_under_trigger_[3]
                                   << " " << dim_not_under_trigger_[4]
                                   << " " << dim_not_under_trigger_[5]
                                   << " " << dim_not_under_trigger_[6] << "]");

        GDEBUG_CONDITION_STREAM(this->verbose.value(), "dimension limits NOT under trigger [Slice Con Phase Rep Set Ave] = [" 
                                   << " " << dim_limit_not_under_trigger_[1] 
                                   << " " << dim_limit_not_under_trigger_[2] 
                                   << " " << dim_limit_not_under_trigger_[3]
                                   << " " << dim_limit_not_under_trigger_[4]
                                   << " " << dim_limit_not_under_trigger_[5]
                                   << " " << dim_limit_not_under_trigger_[6] << "]");

        return GADGET_OK;
    }

    template <typename T, int D>
    void GenericReconAccumulateImageTriggerGadget<T, D>::set_dimensions_under_trigger()
    {
        dim_under_trigger_.resize(num_of_dimensions_, false);
        dim_not_under_trigger_.resize(num_of_dimensions_, false);

        dim_limit_under_trigger_.resize(num_of_dimensions_, 1);
        dim_limit_not_under_trigger_.resize(num_of_dimensions_, 1);

        // ------------------------------------
        // CHA
        // ------------------------------------
        if (cha_trigger_)
        {
            dim_under_trigger_[0] = true;
            dim_limit_under_trigger_[0] = dimensions_[0];
        }
        else
        {
            dim_not_under_trigger_[0] = true;
            dim_limit_not_under_trigger_[0] = dimensions_[0];
        }

        // ------------------------------------
        // SLC
        // ------------------------------------
        if (slc_trigger_)
        {
            dim_under_trigger_[1] = true;
            dim_limit_under_trigger_[1] = dimensions_[1];
        }
        else
        {
            dim_not_under_trigger_[1] = true;
            dim_limit_not_under_trigger_[1] = dimensions_[1];
        }

        // ------------------------------------
        // CON
        // ------------------------------------
        if (con_trigger_)
        {
            dim_under_trigger_[2] = true;
            dim_limit_under_trigger_[2] = dimensions_[2];
        }
        else
        {
            dim_not_under_trigger_[2] = true;
            dim_limit_not_under_trigger_[2] = dimensions_[2];
        }

        // ------------------------------------
        // PHS
        // ------------------------------------
        if (phs_trigger_)
        {
            dim_under_trigger_[3] = true;
            dim_limit_under_trigger_[3] = dimensions_[3];
        }
        else
        {
            dim_not_under_trigger_[3] = true;
            dim_limit_not_under_trigger_[3] = dimensions_[3];
        }

        // ------------------------------------
        // REP
        // ------------------------------------
        if (rep_trigger_)
        {
            dim_under_trigger_[4] = true;
            dim_limit_under_trigger_[4] = dimensions_[4];
        }
        else
        {
            dim_not_under_trigger_[4] = true;
            dim_limit_not_under_trigger_[4] = dimensions_[4];
        }

        // ------------------------------------
        // SET
        // ------------------------------------
        if (set_trigger_)
        {
            dim_under_trigger_[5] = true;
            dim_limit_under_trigger_[5] = dimensions_[5];
        }
        else
        {
            dim_not_under_trigger_[5] = true;
            dim_limit_not_under_trigger_[5] = dimensions_[5];
        }

        // ------------------------------------
        // AVE
        // ------------------------------------
        if (ave_trigger_)
        {
            dim_under_trigger_[6] = true;
            dim_limit_under_trigger_[6] = dimensions_[6];
        }
        else
        {
            dim_not_under_trigger_[6] = true;
            dim_limit_not_under_trigger_[6] = dimensions_[6];
        }

        imageSentBuffer_.create(dim_limit_under_trigger_);
    }

    template <typename T, int D>
    int GenericReconAccumulateImageTriggerGadget<T, D>::process(GadgetContainerMessage<IsmrmrdImageArray>* m1)
    {
        GDEBUG_CONDITION_STREAM(verbose.value(), "GenericReconAccumulateImageTriggerGadget<T, D>::process(...) starts ... ");

        process_called_times_++;

        IsmrmrdImageArray* recon_res_ = m1->getObjectPtr();

        // find the data role
        std::string dataRole = std::string(recon_res_->meta_[0].as_str(GADGETRON_DATA_ROLE));

        size_t encoding = 0;
        if (recon_res_->meta_[0].length("encoding") > 0)
        {
            encoding = (size_t)recon_res_->meta_[0].as_long("encoding", 0);
            GADGET_CHECK_RETURN(encoding < num_encoding_spaces_, GADGET_FAIL);
        }

        size_t ii;
        size_t RO  = recon_res_->data_.get_size(0);
        size_t E1  = recon_res_->data_.get_size(1);
        size_t E2  = recon_res_->data_.get_size(2);
        size_t CHA = recon_res_->data_.get_size(3);
        size_t N   = recon_res_->data_.get_size(4);
        size_t S   = recon_res_->data_.get_size(5);
        size_t SLC = recon_res_->data_.get_size(6);

        image_counter_ += CHA*N*S*SLC;

        GDEBUG_CONDITION_STREAM(this->verbose.value(), "--> receive image array [RO E1 E2 CHA N S SLC] : [" << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC << "]" << " -- " << dataRole);

        if ( (dataRole == GADGETRON_IMAGE_REGULAR) 
            || (dataRole == GADGETRON_IMAGE_RETRO) 
            || (dataRole == GADGETRON_IMAGE_MOCORECON) 
            || (dataRole == GADGETRON_IMAGE_PHASE) 
            || (dataRole == GADGETRON_IMAGE_AIF) )
        {
            // first time call, allocate buffer array
            if(imageBuffer_.get_number_of_elements()==0)
            {
                dimensions_[0] = CHA;
                imageBuffer_.create(dimensions_);
                imageSent_.create(dimensions_);

                size_t nElem = imageBuffer_.get_number_of_elements();
                for (ii = 0; ii<nElem; ii++)
                {
                    imageBuffer_(ii) = ImageType();
                    imageSent_(ii) = false;
                }
            }

            GADGET_CHECK_RETURN(this->store_image(*recon_res_, dimensions_, imageBuffer_)==GADGET_OK, GADGET_FAIL);
            GADGET_CHECK_RETURN(this->trigger(imageBuffer_, imageSent_, false) == GADGET_OK, GADGET_FAIL);
        }

        if ( dataRole == GADGETRON_IMAGE_OTHER )
        {
            if (otherBuffer_.get_number_of_elements() == 0)
            {
                dimensions_[0] = CHA;
                otherBuffer_.create(dimensions_);
                imageSent_.create(dimensions_);

                size_t nElem = otherBuffer_.get_number_of_elements();
                for (ii = 0; ii<nElem; ii++)
                {
                    otherBuffer_(ii) = ImageType();
                    otherSent_(ii) = false;
                }
            }

            GADGET_CHECK_RETURN(this->store_image(*recon_res_, dimensions_, otherBuffer_) == GADGET_OK, GADGET_FAIL);
            GADGET_CHECK_RETURN(this->trigger(otherBuffer_, otherSent_, false) == GADGET_OK, GADGET_FAIL);
        }

        if ( (dataRole==GADGETRON_IMAGE_GFACTOR) 
            || (dataRole==GADGETRON_IMAGE_SNR_MAP) 
            || (dataRole==GADGETRON_IMAGE_STD_MAP) 
            || (dataRole==GADGETRON_IMAGE_WRAPAROUNDMAP)
            || (dataRole == GADGETRON_IMAGE_RECON_FIGURE)
            || pass_image_immediate_ )
        {
            std::vector<size_t> dimIm(D,0);
            dimIm[0] = RO;
            if (D>1) dimIm[1] = E1;
            if (D>2) dimIm[2] = E2;

            std::vector<size_t> dim(num_of_dimensions_, 1);

            // pass the image to the next gadget
            size_t slc, n, s, cha;
            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        for (cha = 0; cha < CHA; cha++)
                        {
                            Gadgetron::GadgetContainerMessage<ImageBufferType>* cm1 = new Gadgetron::GadgetContainerMessage<ImageBufferType>();
                            ImageBufferType& imgBuf = *(cm1->getObjectPtr());
                            imgBuf.create(dim);

                            // set image content
                            imgBuf(0).create(dimIm);
                            ValueType* pIm = imgBuf(0).begin();
                            size_t numPixel = imgBuf(0).get_number_of_elements();

                            std::complex<float>* pData = &(recon_res_->data_(0, 0, 0, cha, n, s, slc));

                            size_t ii;
                            for (ii = 0; ii < numPixel; ii++)
                            {
                                pIm[ii] = (ValueType)(pData[ii]);
                            }

                            // set image attrib
                            imgBuf(0).header_ = recon_res_->headers_(n, s, slc);
                            imgBuf(0).attrib_ = recon_res_->meta_[n + s*N + slc*N*S];

                            imgBuf(0).set_pixel_size(0, imgBuf(0).header_.field_of_view[0] / RO);
                            if (D>1) imgBuf(0).set_pixel_size(1, imgBuf(0).header_.field_of_view[1] / E1);
                            if (D>2) imgBuf(0).set_pixel_size(2, imgBuf(0).header_.field_of_view[2] / E2);

                            imgBuf(0).set_image_position(0, imgBuf(0).header_.position[0]);
                            imgBuf(0).set_image_position(1, imgBuf(0).header_.position[1]);
                            imgBuf(0).set_image_position(2, imgBuf(0).header_.position[2]);

                            imgBuf(0).set_image_orientation(0, 0, imgBuf(0).header_.read_dir[0]);
                            imgBuf(0).set_image_orientation(0, 1, imgBuf(0).header_.read_dir[1]);
                            imgBuf(0).set_image_orientation(0, 2, imgBuf(0).header_.read_dir[2]);

                            imgBuf(0).set_image_orientation(1, 0, imgBuf(0).header_.phase_dir[0]);
                            imgBuf(0).set_image_orientation(1, 1, imgBuf(0).header_.phase_dir[1]);
                            imgBuf(0).set_image_orientation(1, 2, imgBuf(0).header_.phase_dir[2]);

                            imgBuf(0).set_image_orientation(2, 0, imgBuf(0).header_.slice_dir[0]);
                            imgBuf(0).set_image_orientation(2, 1, imgBuf(0).header_.slice_dir[1]);
                            imgBuf(0).set_image_orientation(2, 2, imgBuf(0).header_.slice_dir[2]);

                            imgBuf(0).attrib_.set(GADGETRON_PASS_IMMEDIATE, (long)1);

                            if (this->next()->putq(cm1) < 0)
                            {
                                m1->release();
                                cm1->release();
                                return GADGET_FAIL;
                            }
                        }
                    }
                }
            }
        }

        m1->release();
        return GADGET_OK;
    }

    template <typename T, int D>
    int GenericReconAccumulateImageTriggerGadget<T, D>::trigger(ImageBufferType& buf, ImageSentFlagBufferType& sentFlagBuf, bool inClose)
    {
        try
        {
            // scan the buffered images, if the trigger dimensions are complete, sent out this package
            if(buf.get_number_of_elements()==0 || !buf.dimensions_equal(this->dimensions_))
            {
                return GADGET_OK;
            }

            // not under trigger
            size_t cha, slc, con, phs, rep, set, ave;

            // under trigger
            size_t cha_t, slc_t, con_t, phs_t, rep_t, set_t, ave_t;

            std::vector<size_t> image_ind(num_of_dimensions_, 0);
            std::vector<size_t> image_sent_ind(num_of_dimensions_, 0);

            size_t numOfElem = imageSentBuffer_.get_number_of_elements();
            size_t ii;
            for ( ii=0; ii<numOfElem; ii++ ) { imageSentBuffer_(ii) = ImageType(); }

            for ( ave=0; ave<dim_limit_not_under_trigger_[6]; ave++ )
            {
                if ( dim_not_under_trigger_[6] ) image_ind[6] = ave;
                // -------------------
                for ( set=0; set<dim_limit_not_under_trigger_[5]; set++ )
                {
                    if ( dim_not_under_trigger_[5] ) image_ind[5] = set;
                    // -------------------
                    for ( rep=0; rep<dim_limit_not_under_trigger_[4]; rep++ )
                    {
                        if ( dim_not_under_trigger_[4] ) image_ind[4] = rep;
                        // -------------------
                        for ( phs=0; phs<dim_limit_not_under_trigger_[3]; phs++ )
                        {
                            if ( dim_not_under_trigger_[3] ) image_ind[3] = phs;
                            // -------------------
                            for ( con=0; con<dim_limit_not_under_trigger_[2]; con++ )
                            {
                                if ( dim_not_under_trigger_[2] ) image_ind[2] = con;
                                // -------------------
                                for ( slc=0; slc<dim_limit_not_under_trigger_[1]; slc++ )
                                {
                                    if ( dim_not_under_trigger_[1] ) image_ind[1] = slc;
                                    // -------------------
                                    for ( cha=0; cha<dim_limit_not_under_trigger_[0]; cha++ )
                                    {
                                        if ( dim_not_under_trigger_[0] ) image_ind[0] = cha;
                                        // -------------------

                                        // loop over under triggered dimensions and check whether every images are there
                                        bool needTrigger = true;
                                        if ( inClose )
                                        {
                                            needTrigger = false;
                                        }

                                        {
                                            for ( ii=0; ii<numOfElem; ii++ ) { imageSentBuffer_(ii) = ImageType(); }

                                            // =================================================

                                            for ( ave_t=0; ave_t<dim_limit_under_trigger_[6]; ave_t++ )
                                            {
                                                if ( dim_under_trigger_[6] ) image_ind[6] = ave_t;
                                                image_sent_ind[6] = ave_t;
                                                // -------------------
                                                for ( set_t=0; set_t<dim_limit_under_trigger_[5]; set_t++ )
                                                {
                                                    if ( dim_under_trigger_[5] ) image_ind[5] = set_t;
                                                    image_sent_ind[5] = set_t;
                                                    // -------------------
                                                    for ( rep_t=0; rep_t<dim_limit_under_trigger_[4]; rep_t++ )
                                                    {
                                                        if ( dim_under_trigger_[4] ) image_ind[4] = rep_t;
                                                        image_sent_ind[4] = rep_t;
                                                        // -------------------
                                                        for ( phs_t=0; phs_t<dim_limit_under_trigger_[3]; phs_t++ )
                                                        {
                                                            if ( dim_under_trigger_[3] ) image_ind[3] = phs_t;
                                                            image_sent_ind[3] = phs_t;
                                                            // -------------------
                                                            for ( con_t=0; con_t<dim_limit_under_trigger_[2]; con_t++ )
                                                            {
                                                                if ( dim_under_trigger_[2] ) image_ind[2] = con_t;
                                                                image_sent_ind[2] = con_t;
                                                                // -------------------
                                                                for ( slc_t=0; slc_t<dim_limit_under_trigger_[1]; slc_t++ )
                                                                {
                                                                    if ( dim_under_trigger_[1] ) image_ind[1] = slc_t;
                                                                    image_sent_ind[1] = slc_t;
                                                                    // -------------------
                                                                    for ( cha_t=0; cha_t<dim_limit_under_trigger_[0]; cha_t++ )
                                                                    {
                                                                        if ( dim_under_trigger_[0] ) image_ind[0] = cha_t;
                                                                        image_sent_ind[0] = cha_t;
                                                                        // -------------------

                                                                        ImageType& img = buf(image_ind);
                                                                        bool sentFlag = sentFlagBuf(image_ind);

                                                                        if ( inClose )
                                                                        {
                                                                            // if in close call, send out all unsent images
                                                                            if ( img.get_number_of_elements()>0 && !sentFlag )
                                                                            {
                                                                                imageSentBuffer_(image_sent_ind) = std::move(img);
                                                                                buf(image_ind) = ImageType();
                                                                                needTrigger = true;
                                                                            }
                                                                        }
                                                                        else
                                                                        {
                                                                            if (img.get_number_of_elements()>0 && !sentFlag )
                                                                            {
                                                                                imageSentBuffer_(image_sent_ind) = img;
                                                                            }
                                                                            else
                                                                            {
                                                                                needTrigger = false; // if all images for current under-trigger dimensions are filled, trigger
                                                                                break;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }

                                            if ( needTrigger )
                                            {
                                                // if a image has been sent, not sent again
                                                for ( ave_t=0; ave_t<dim_limit_under_trigger_[6]; ave_t++ )
                                                {
                                                    if ( dim_under_trigger_[6] ) image_ind[6] = ave_t;
                                                    // -------------------
                                                    for ( set_t=0; set_t<dim_limit_under_trigger_[5]; set_t++ )
                                                    {
                                                        if ( dim_under_trigger_[5] ) image_ind[5] = set_t;
                                                        // -------------------
                                                        for ( rep_t=0; rep_t<dim_limit_under_trigger_[4]; rep_t++ )
                                                        {
                                                            if ( dim_under_trigger_[4] ) image_ind[4] = rep_t;
                                                            // -------------------
                                                            for ( phs_t=0; phs_t<dim_limit_under_trigger_[3]; phs_t++ )
                                                            {
                                                                if ( dim_under_trigger_[3] ) image_ind[3] = phs_t;
                                                                // -------------------
                                                                for ( con_t=0; con_t<dim_limit_under_trigger_[2]; con_t++ )
                                                                {
                                                                    if ( dim_under_trigger_[2] ) image_ind[2] = con_t;
                                                                    // -------------------
                                                                    for ( slc_t=0; slc_t<dim_limit_under_trigger_[1]; slc_t++ )
                                                                    {
                                                                        if ( dim_under_trigger_[1] ) image_ind[1] = slc_t;
                                                                        // -------------------
                                                                        for ( cha_t=0; cha_t<dim_limit_under_trigger_[0]; cha_t++ )
                                                                        {
                                                                            if ( dim_under_trigger_[0] ) image_ind[0] = cha_t;
                                                                            // -------------------

                                                                            bool sentFlag = sentFlagBuf(image_ind);
                                                                            if ( sentFlag )
                                                                            {
                                                                                imageSentBuffer_(cha_t, slc_t, con_t, phs_t, rep_t, set_t, ave_t) = ImageType();
                                                                            }
                                                                            else
                                                                            {
                                                                                sentFlagBuf(image_ind) = true;
                                                                            }

                                                                            buf(image_ind) = ImageType();
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }

                                                GDEBUG_CONDITION_STREAM(verbose.value(), "--> Accumulator image trigger for [CHA SLC CON PHS REP SET AVE] : ["
                                                                                                                            << image_ind[0] << " " 
                                                                                                                            << image_ind[1] << " " 
                                                                                                                            << image_ind[2] << " " 
                                                                                                                            << image_ind[3] << " " 
                                                                                                                            << image_ind[4] << " " 
                                                                                                                            << image_ind[5] << " " 
                                                                                                                            << image_ind[6] << "]" );

                                                Gadgetron::GadgetContainerMessage<ImageBufferType>* cm1 = new Gadgetron::GadgetContainerMessage<ImageBufferType>();
                                                ImageBufferType& imgBuf = *(cm1->getObjectPtr());
                                                imgBuf = imageSentBuffer_;

                                                if (!this->wave_form_buffer_.empty())
                                                {
                                                    Gadgetron::GadgetContainerMessage<std::vector<ISMRMRD::Waveform>>* cm2 = new Gadgetron::GadgetContainerMessage<std::vector<ISMRMRD::Waveform>>();

                                                    *(cm2->getObjectPtr()) = this->wave_form_buffer_;
                                                    this->wave_form_buffer_.clear();

                                                    cm1->cont(cm2);
                                                }

                                                if (this->next()->putq(cm1) < 0) 
                                                {
                                                    cm1->release();
                                                    return GADGET_FAIL;
                                                }
                                            }
                                            else
                                            {
                                                for ( ii=0; ii<numOfElem; ii++ )
                                                {
                                                    imageSentBuffer_(ii) = ImageType();
                                                }
                                            }

                                            // =================================================
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
            GERROR_STREAM("Error happens in GenericReconAccumulateImageTriggerGadget<T, D>::trigger(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    template <typename T, int D>
    int GenericReconAccumulateImageTriggerGadget<T, D>::store_image(const IsmrmrdImageArray& img, const std::vector<size_t>& buf_dimension, ImageBufferType& buf)
    {
        try
        {
            size_t RO  = img.data_.get_size(0);
            size_t E1  = img.data_.get_size(1);
            size_t E2  = img.data_.get_size(2);
            size_t CHA = img.data_.get_size(3);
            size_t N   = img.data_.get_size(4);
            size_t S   = img.data_.get_size(5);
            size_t SLC = img.data_.get_size(6);

            std::vector<size_t> dimIm(D, 0);
            dimIm[0] = RO;
            if (D>1) dimIm[1] = E1;
            if (D>2) dimIm[2] = E2;

            size_t slc, s, n, cha;

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        for (cha = 0; cha < CHA; cha++)
                        {
                            size_t slice = img.headers_(n, s, slc).slice;
                            size_t con   = img.headers_(n, s, slc).contrast;
                            size_t phs   = img.headers_(n, s, slc).phase;
                            size_t rep   = img.headers_(n, s, slc).repetition;
                            size_t set   = img.headers_(n, s, slc).set;
                            size_t ave   = img.headers_(n, s, slc).average;

                            if ((cha >= buf_dimension[0])
                                || (slice >= buf_dimension[1])
                                || (con >= buf_dimension[2])
                                || (phs >= buf_dimension[3])
                                || (rep >= buf_dimension[4])
                                || (set >= buf_dimension[5])
                                || (ave >= buf_dimension[6]))
                            {
                                GERROR_STREAM("Incoming image is out of boundary [cha slc con phs rep set ave] : [" << cha << " " << slice << " " << con << " " << phs << " " << rep << " " << set << " " << ave << "]");
                                return GADGET_FAIL;
                            }

                            // create image
                            ImageType storedImage;
                            storedImage.create(dimIm);

                            ValueType* pIm = storedImage.begin();
                            size_t numPixel = storedImage.get_number_of_elements();

                            const std::complex<float>* pData = &(img.data_(0, 0, 0, cha, n, s, slc));
                            for (size_t ii = 0; ii < numPixel; ii++)
                            {
                                pIm[ii] = (ValueType)(pData[ii]);
                            }

                            storedImage.attrib_ = img.meta_[n + s*N + slc*N*S];
                            storedImage.header_ = img.headers_(n, s, slc);

                            storedImage.set_pixel_size(0, img.headers_(n, s, slc).field_of_view[0] / RO);
                            if (D>1) storedImage.set_pixel_size(1, img.headers_(n, s, slc).field_of_view[1] / E1);
                            if (D>2) storedImage.set_pixel_size(2, img.headers_(n, s, slc).field_of_view[2] / E2);

                            storedImage.set_image_position(0, img.headers_(n, s, slc).position[0]);
                            storedImage.set_image_position(1, img.headers_(n, s, slc).position[1]);
                            storedImage.set_image_position(2, img.headers_(n, s, slc).position[2]);

                            storedImage.set_image_orientation(0, 0, img.headers_(n, s, slc).read_dir[0]);
                            storedImage.set_image_orientation(0, 1, img.headers_(n, s, slc).read_dir[1]);
                            storedImage.set_image_orientation(0, 2, img.headers_(n, s, slc).read_dir[2]);

                            storedImage.set_image_orientation(1, 0, img.headers_(n, s, slc).phase_dir[0]);
                            storedImage.set_image_orientation(1, 1, img.headers_(n, s, slc).phase_dir[1]);
                            storedImage.set_image_orientation(1, 2, img.headers_(n, s, slc).phase_dir[2]);

                            storedImage.set_image_orientation(2, 0, img.headers_(n, s, slc).slice_dir[0]);
                            storedImage.set_image_orientation(2, 1, img.headers_(n, s, slc).slice_dir[1]);
                            storedImage.set_image_orientation(2, 2, img.headers_(n, s, slc).slice_dir[2]);

                            storedImage.attrib_.set(GADGETRON_PASS_IMMEDIATE, (long)0);

                            buf(cha, slice, con, phs, rep, set, ave) = storedImage;
                        }
                    }
                }
            }

            // store the waveform
            if (img.waveform_.has_value())
            {
                for (size_t ii = 0; ii < img.waveform_.value().size(); ii++)
                {
                    wave_form_buffer_.push_back((*img.waveform_)[ii]);
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happens in GenericReconAccumulateImageTriggerGadget<T, D>::storeImage(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    template <typename T, int D> 
    int GenericReconAccumulateImageTriggerGadget<T, D>::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GenericReconAccumulateImageTriggerGadget - close(flags) : " << flags);

        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

        if ( flags!=0 && !triggered_in_close_ )
        {
            triggered_in_close_ = true;

            GDEBUG_CONDITION_STREAM(true, "GenericReconAccumulateImageTriggerGadget - trigger in close(flags) ... ");

            GADGET_CHECK_RETURN(this->trigger(imageBuffer_, imageSent_, true) == GADGET_OK, GADGET_FAIL);
            GADGET_CHECK_RETURN(this->trigger(otherBuffer_, otherSent_, true) == GADGET_OK, GADGET_FAIL);
        }

        return GADGET_OK;
    }

    GenericReconAccumulateImage2DTriggerGadget::GenericReconAccumulateImage2DTriggerGadget() : BaseClass()
    {
    }

    GenericReconAccumulateImage2DTriggerGadget::~GenericReconAccumulateImage2DTriggerGadget()
    {
    }

    GenericReconAccumulateImage3DTriggerGadget::GenericReconAccumulateImage3DTriggerGadget() : BaseClass()
    {
    }

    GenericReconAccumulateImage3DTriggerGadget::~GenericReconAccumulateImage3DTriggerGadget()
    {
    }

    GADGET_FACTORY_DECLARE(GenericReconAccumulateImage2DTriggerGadget)
    GADGET_FACTORY_DECLARE(GenericReconAccumulateImage3DTriggerGadget)
}
